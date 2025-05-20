#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <chrono>
#include <unordered_map>
#include <memory>
#include <limits>
#include <filesystem>
#include <unistd.h>

#include "compressorData.hpp"
#include "evtk.hpp"
#include "scalarfield.hpp"
#include "scalarfieldheap.hpp"
#include "unionfind.hpp"
#include "multitracker.hpp"
#include "unrolledloops.hpp"
#include "criticalpointstorage.hpp"
#include "criticalpointtree.hpp"
#include "mergeTreeHelper.hpp"
#include "unorderedunionfind.hpp"
#include "region.hpp"
#include "huffman.hpp"
#include "getCTFull.hpp"
#include "cubicsplineinterpolation.hpp"
#include "results.hpp"
#include "write_utils.hpp"
#include "python_error.hpp"

using namespace std;

template <typename T> 
Results* compressIterativeMT(string filename, string arrayName, double epsilon_relative, double xi_relative, string outputFilename = "compressed", string outputFolder = ".", string baseCompressor = "SZ3", double compressorParameter = 1, bool verbose = false, bool logQuantizeMode = true, int initialPrecision = 0, string baseCompressorFolder = "../../base_compressors", int size_x = -1, int size_y = -1, int size_z = -1 ){
    Results* results = new Results();

    auto startTime = std::chrono::high_resolution_clock::now();

    ScalarField<T> sf;
    if( size_z == -1 ){
        sf.loadFromVTK(filename, arrayName);
    } else {
        sf.loadFromDat(filename, size_x, size_y, size_z, true);
    }

    int numPoints = sf.size();

    T epsilon = epsilon_relative * sf.getDataRange();
    T xi = xi_relative * sf.getDataRange();

    string systemSuffix = "";
    if( !verbose ){
        systemSuffix = " > /dev/null 2> /dev/null";
    }

    // initial compression
    if( baseCompressor == "SZ3" ){

        string typeFlag;
        if( is_same<T, double>::value ){
            typeFlag = "-d";
        } else if( is_same<T, float>::value ){
            typeFlag = "-f";
        }

        sf.saveToDat(outputFolder + "/rawData.dat");
        int result = system( (baseCompressorFolder + "/sz3 " + typeFlag + " -i " + outputFolder + "/rawData.dat -z " + outputFolder + "/rawData.cmp -3 " + to_string(sf.size_x) + " " + to_string(sf.size_y) + " " + to_string(sf.size_z) + " -M ABS " + to_string(compressorParameter*xi) + systemSuffix).c_str() );
        result = system( (baseCompressorFolder + "/sz3 " + typeFlag + " -z " + outputFolder + "/rawData.cmp -o " + outputFolder + "/intermediateData.dat -3 " + to_string(sf.size_x) + " " + to_string(sf.size_y) + " " + to_string(sf.size_z) + " -M ABS " + to_string(compressorParameter*xi) + systemSuffix).c_str() );

    } else if( baseCompressor == "ZFP" ){

        string typeFlag;
        if( is_same<T, double>::value ){
            typeFlag = "-d";
        } else if( is_same<T, float>::value ){
            typeFlag = "-f";
        }

        sf.saveToDat(outputFolder + "/rawData.dat");
        int result = system( (baseCompressorFolder + "/zfp " + typeFlag + " -i " + outputFolder + "/rawData.dat -z " + outputFolder + "/rawData.cmp -3 " + to_string(sf.size_x) + " " + to_string(sf.size_y) + " " + to_string(sf.size_z) + " -a " + to_string(compressorParameter*xi) + systemSuffix).c_str());
        result = system( (baseCompressorFolder + "/zfp " + typeFlag + " -z " + outputFolder + "/rawData.cmp -o " + outputFolder + "/intermediateData.dat -3 " + to_string(sf.size_x) + " " + to_string(sf.size_y) + " " + to_string(sf.size_z) + " -a " + to_string(compressorParameter*xi) + systemSuffix).c_str() );

    } else if( baseCompressor == "TTHRESH" ){

        string typeFlag;
        if( is_same<T, double>::value ){
            typeFlag = "double";
        } else if( is_same<T, float>::value ){
            typeFlag = "float";
        }

        sf.saveToDat(outputFolder + "/rawData.dat");
        int result = system( (baseCompressorFolder + "/tthresh -i " + outputFolder + "/rawData.dat -t " + typeFlag + " -s " + to_string(sf.size_x) + " " + to_string(sf.size_y) + " " + to_string(sf.size_z) + " -r " + to_string(compressorParameter*xi) + " -c " + outputFolder + "/rawData.cmp -o " + outputFolder + "/intermediateData.dat" + systemSuffix).c_str());
    } else if( baseCompressor == "CSI" ){
        int ratio = lround(compressorParameter);
        ScalarField<float>* downsampled = sf.downSample();
        csiCompressScalarField<float>(*downsampled, outputFolder + "/rawData.cmp", ratio);
        delete downsampled;
        ScalarField<float>* reconstructed32 = csiDecompress<float>(outputFolder + "/rawData.cmp", sf.size_x, sf.size_y, sf.size_z, ratio);

        if( is_same<T,float>::value ){
            reconstructed32->saveToDat(outputFolder + "/intermediateData.dat");
            delete reconstructed32;
        } else {
            ScalarField<double>* reconstructed64 = reconstructed32->upSample();
            reconstructed64->saveToDat(outputFolder + "/intermediateData.dat");
            delete reconstructed32;
            delete reconstructed64;

        }

    } else if( baseCompressor == "Neurcomp" ){

        int clusterBits = lround(compressorParameter);

        int result = system(("python3 " + baseCompressorFolder + "/neurcomp/net_compress.py --net " + outputFolder + "/thenet.pth --config " + outputFolder + "/thenet.json --compressed " + outputFolder + "/rawData.cmp --cluster_bits " + to_string(clusterBits) + systemSuffix).c_str());
        pythonError(result, verbose);

        result = system( ("python3 " + baseCompressorFolder + "/neurcomp/net_decompress.py --volume " + outputFolder + "/dummy.npy --compressed " + outputFolder + "/rawData.cmp --recon " + outputFolder + "/recon" + systemSuffix).c_str() );
        pythonError(result, verbose);

        ScalarField<T> reconVTI;
        reconVTI.loadFromVTK(outputFolder + "/recon.vti", "sf");
        reconVTI.saveToDat(outputFolder + "/intermediateData.dat");
        result = system(("rm " + outputFolder + "/recon.vti" + systemSuffix).c_str());
    }

    auto baseCompressionSplit = std::chrono::high_resolution_clock::now();

    ScalarField<T> sfIntermediate;
    sfIntermediate.loadFromDat( outputFolder + "/intermediateData.dat", sf.size_x, sf.size_y, sf.size_z );

    ContourTreeInfo<T> contourTreeInfo = getContourTree(&sf, epsilon_relative, verbose);

    // delete unneeded things.
    for( auto i : *(contourTreeInfo.contourTreeNodes) ){
        contourTreeInfo.contourTree->deleteNode(i);
    }
    delete contourTreeInfo.contourTree;
    delete contourTreeInfo.contourTreeNodes;
    contourTreeInfo.contourTree = NULL;
    contourTreeInfo.contourTreeNodes = NULL;

    MergeTreeInfo<T>* joinTreeInfo = contourTreeInfo.joinTreeInfo;
    MergeTreeInfo<T>* splitTreeInfo = contourTreeInfo.splitTreeInfo;

    auto contourTreeSplit = std::chrono::high_resolution_clock::now();

    // pair<MergeTreeInfo<T>*,MergeTreeInfo<T>*> mti = getJoinAndSplitTrees<T>(&sf, epsilon, verbose);
    // MergeTreeInfo<T>* joinTreeInfo = mti.first;
    // MergeTreeInfo<T>* splitTreeInfo = mti.second;

    // create empty scalar field that will store the final values.
    ScalarField<T> sfFinal = ScalarField<T>();
    sfFinal.createEmptySF(sf.size_x, sf.size_y, sf.size_z);

    // upper bound, lower bound, quantization, precision storage are abstracted
    unique_ptr<CompressorData<T>> cd(new CompressorData<T>(&sf, &sfIntermediate, &sfFinal, xi, logQuantizeMode, initialPrecision));

    // find the segmentation of each points. Set initial error bounds based on segmentation and store the segmentation.

    // maps that will store lists of each point needed for segmentation
    unordered_map<int, unique_ptr<vector<int>>> jtBaseSegmentationLists;
    unordered_map<int, unique_ptr<vector<int>>> stBaseSegmentationLists;

    // comparisons needed for the following step
    ScalarFieldMore<T> groundIndexIsMore;
    ScalarFieldLess<T> groundIndexIsLess;

    groundIndexIsMore.setSF(&sf);
    groundIndexIsLess.setSF(&sf);

    // aliased variables to make things nicer
    vector<int>& jtUSSC = *(joinTreeInfo->USSC);
    vector<int>& jtUSCI = *(joinTreeInfo->USCI);
    vector<int>& stUSSC = *(splitTreeInfo->USSC);
    vector<int>& stUSCI = *(splitTreeInfo->USCI);

    unordered_map<int,int>& jtVTTS = *(joinTreeInfo->visitedTrackerIdToSegmentation);
    unordered_map<int,int>& stVTTS = *(splitTreeInfo->visitedTrackerIdToSegmentation);

    // used for tracking branches so that we can adjust segmentations on the fly to make sure that they are valid.
    // first entry of the map is the direction to the smallest child. Second is the smallest child itself.
    unordered_map<int,pair<int,int>>* jtPathToSmallestVertexPtr = new unordered_map<int,pair<int,int>>();
    unordered_map<int,pair<int,int>>* stPathToLargestVertexPtr = new unordered_map<int,pair<int,int>>();

    unordered_map<int,pair<int,int>>& jtPathToSmallestVertex = *jtPathToSmallestVertexPtr;
    unordered_map<int,pair<int,int>>& stPathToLargestVertex = *stPathToLargestVertexPtr;

    unordered_map<int,int>* numChildrenChecked = new unordered_map<int,int>();
    vector<int>* verticesToExpand = new vector<int>();

    for( auto node : *(joinTreeInfo->vertices) ){
        if( node->criticalType == 0 ){
            jtPathToSmallestVertex[node->globalIndex] = { node->globalIndex, node->globalIndex };
            verticesToExpand->push_back(node->globalIndex);
        }
    }

    for( int i = 0; i < verticesToExpand->size(); ++i ){
        int nodeIndex = (*verticesToExpand)[i];
        int parentIndex = joinTreeInfo->tree->getParentGlobalIndex(nodeIndex);
        if( parentIndex != -1 ){
            if( !jtPathToSmallestVertex.count(parentIndex) || groundIndexIsMore(jtPathToSmallestVertex[parentIndex].second, jtPathToSmallestVertex[nodeIndex].second) ){
                jtPathToSmallestVertex[parentIndex] = { nodeIndex, jtPathToSmallestVertex[nodeIndex].second };
            }
            if( !numChildrenChecked->count(parentIndex) ){
                (*numChildrenChecked)[parentIndex] = 1;
            } else {
                ++((*numChildrenChecked)[parentIndex]);
            }

            if( (*numChildrenChecked)[parentIndex] == joinTreeInfo->tree->getNode(parentIndex)->children.size() ){
                verticesToExpand->push_back(parentIndex);
            }
        }
    }

    // repeat but for split tree

    delete numChildrenChecked;
    
    numChildrenChecked = new unordered_map<int,int>();    
    verticesToExpand->clear();

    for( auto node : *(splitTreeInfo->vertices) ){
        if( node->criticalType == 3 ){
            stPathToLargestVertex[node->globalIndex] = { node->globalIndex, node->globalIndex };
            verticesToExpand->push_back(node->globalIndex);
        }
    }

    for( int i = 0; i < verticesToExpand->size(); ++i ){
        int nodeIndex = (*verticesToExpand)[i];
        int parentIndex = splitTreeInfo->tree->getParentGlobalIndex(nodeIndex);
        if( parentIndex != -1 ){
            if( !stPathToLargestVertex.count(parentIndex) || groundIndexIsLess(stPathToLargestVertex[parentIndex].second, stPathToLargestVertex[nodeIndex].second) ){
                stPathToLargestVertex[parentIndex] = { nodeIndex, stPathToLargestVertex[nodeIndex].second };
            }
            if( !numChildrenChecked->count(parentIndex) ){
                (*numChildrenChecked)[parentIndex] = 1;
            } else {
                ++((*numChildrenChecked)[parentIndex]);
            }

            if( (*numChildrenChecked)[parentIndex] == splitTreeInfo->tree->getNode(parentIndex)->children.size() ){
                verticesToExpand->push_back(parentIndex);
            }
        }
    }

    delete numChildrenChecked;
    delete verticesToExpand;

    int badEB = 0;

    // iterate through all points. Calculate bounds based on segmentation and quantize each point.
    for(int i = 0; i < numPoints; ++i){

        if( joinTreeInfo->tree->hasNode(i) || splitTreeInfo->tree->hasNode(i) ){
            // store critical points losslessly

            if( joinTreeInfo->tree->hasNode(i) && joinTreeInfo->tree->getNode(i)->criticalType == 0 ){
                if( !jtBaseSegmentationLists.count( i ) ){
                    jtBaseSegmentationLists[i] = unique_ptr<vector<int>>(new vector<int>());
                }
                jtBaseSegmentationLists[i]->push_back(i);
            } else if( splitTreeInfo->tree->hasNode(i) && splitTreeInfo->tree->getNode(i)->criticalType == 3 ){
                if( !stBaseSegmentationLists.count(i) ){
                    stBaseSegmentationLists[i] = unique_ptr<vector<int>>(new vector<int>());
                }
                stBaseSegmentationLists[i]->push_back(i);
            }

            cd->initializeNextPointLossless();

        } else {

            vector<int> jtSegmentations;
            vector<int> stSegmentations;

            // calculate error bound for regular points

            // first, identify the join tree and split tree segmentation

            // bounds from join tree
            int jtSegmentation = joinTreeInfo->visitedTracker->entries[i];
            if( jtSegmentation == -1 ){

                // in this case, the join tree segmentation comes from inside of saddles.
                jtSegmentations.push_back( jtVTTS[ jtUSCI[upper_bound( jtUSSC.begin(), jtUSSC.end(), i, groundIndexIsLess ) - jtUSSC.begin() - 1] ] );
            } else if( jtSegmentation == -2 ){

                for( auto subSegmentation : joinTreeInfo->visitedTracker->spillOverMarks[ i ] ){
                    jtSegmentations.push_back( jtVTTS[subSegmentation] );
                }

            } else { 
                jtSegmentations.push_back( jtVTTS[ jtSegmentation ] );
            }

            // bounds from split tree
            int stSegmentation = splitTreeInfo->visitedTracker->entries[i];
            if( stSegmentation == -1 ){
                // in this case, the join tree segmentation comes from inside of saddles.
                stSegmentations.push_back(stVTTS[ stUSCI[upper_bound( stUSSC.begin(), stUSSC.end(), i, groundIndexIsMore ) - stUSSC.begin() - 1] ]);
            } else if( stSegmentation == -2 ){
                for( auto subSegmentation : splitTreeInfo->visitedTracker->spillOverMarks[ i ] ){
                    stSegmentations.push_back(stVTTS[subSegmentation]);
                }
            } else { 
                stSegmentations.push_back(stVTTS[stSegmentation]);
            }

            // then add this point to all relevant segmentation lists
            for( auto jtSegmentationElt : jtSegmentations ){
                if( !jtBaseSegmentationLists.count(jtSegmentationElt) ){
                    jtBaseSegmentationLists[jtSegmentationElt] = unique_ptr<vector<int>>(new vector<int>());
                }
                jtBaseSegmentationLists[jtSegmentationElt]->push_back(i);
            }

            for( auto stSegmentationElt : stSegmentations ){
                if( !stBaseSegmentationLists.count(stSegmentationElt) ){
                    stBaseSegmentationLists[stSegmentationElt] = unique_ptr<vector<int>>(new vector<int>());
                }
                stBaseSegmentationLists[stSegmentationElt]->push_back(i);
            }

            // If a point's value is not within the range of its stated segmentation, which can be the result
            // of merging in persistent simplification, move down the branch until the segmentation is valid.
            for( int j = 0; j < jtSegmentations.size(); ++j ){
                int segmentation = jtSegmentations[j];
                while( groundIndexIsMore(segmentation, i) ){
                    segmentation = jtPathToSmallestVertex[segmentation].first;
                }
                jtSegmentations[j] = segmentation;
            }

            for( int j = 0; j < stSegmentations.size(); ++j ){
                int segmentation = stSegmentations[j];
                while( groundIndexIsLess(segmentation, i) ){
                    segmentation = stPathToLargestVertex[segmentation].first;
                }
                stSegmentations[j] = segmentation;
            }

            // set bounds from combined bounds
            T lowerBound = sf.getElement(i)-xi;
            T upperBound = sf.getElement(i)+xi;

            for( auto jtSegmentationElt : jtSegmentations ){
                for( auto stSegmentationElt : stSegmentations ){
                    pair<T,T> segmentationBounds = (*(contourTreeInfo.segmentBounds))[jtSegmentationElt][stSegmentationElt];

                    lowerBound = max( lowerBound, segmentationBounds.first );
                    upperBound = min( upperBound, segmentationBounds.second );
                }
            }

            cd->initializeNextPointFromBounds(lowerBound, upperBound);
        }

    }

    if( verbose ){
        cout << "set initial error bounds" << endl;
    }

    auto errorBoundSplit = std::chrono::high_resolution_clock::now();

    // delete variables that are only needed for the above and take up a lot of space
    delete joinTreeInfo->USCI;
    delete joinTreeInfo->USSC;
    delete joinTreeInfo->visitedTracker;
    delete joinTreeInfo->visitedTrackerIdToSegmentation;
    delete splitTreeInfo->USCI;
    delete splitTreeInfo->USSC;
    delete splitTreeInfo->visitedTracker;
    delete splitTreeInfo->visitedTrackerIdToSegmentation;
    delete contourTreeInfo.segmentBounds;

    joinTreeInfo->USCI = NULL;
    joinTreeInfo->USSC = NULL;
    joinTreeInfo->visitedTracker = NULL;
    joinTreeInfo->visitedTrackerIdToSegmentation = NULL;
    splitTreeInfo->USCI = NULL;
    splitTreeInfo->USSC = NULL;
    splitTreeInfo->visitedTracker = NULL;
    splitTreeInfo->visitedTrackerIdToSegmentation = NULL;
    contourTreeInfo.segmentBounds = NULL;

    bool keepGoing = true;
    int step = 0;
    while( keepGoing ){
        // find false cases

        pair<MergeTreeInfo<T>*,MergeTreeInfo<T>*> reconMti = getJoinAndSplitTrees<T>(cd->reconstructedData, epsilon, verbose);
        MergeTreeInfo<T>* joinTreeInfo2 = reconMti.first;
        MergeTreeInfo<T>* splitTreeInfo2 = reconMti.second;

        // check for fp and fn
        vector<int> falsePositivesJT;
        vector<int> falseNegativesJT;
        vector<int> falsePositivesST;
        vector<int> falseNegativesST;

        // jt false cases
        for( auto node : *(joinTreeInfo->vertices) ){
            if( node->criticalType == 0 ){
                int i = node->globalIndex;
                if ( !joinTreeInfo2->tree->hasNode(i) || joinTreeInfo2->tree->getParentGlobalIndex(i) != joinTreeInfo->tree->getParentGlobalIndex(i)
                     || joinTreeInfo2->tree->getNode(i)->criticalType != node->criticalType ){
                    falseNegativesJT.push_back(i);
                }
            }

        }

        for( auto node : *(joinTreeInfo2->vertices) ){
            if( node->criticalType == 0 ){
                int i = node->globalIndex;
                if ( !joinTreeInfo->tree->hasNode(i) || joinTreeInfo->tree->getParentGlobalIndex(i) != node->parent->globalIndex ) {
                    falsePositivesJT.push_back(i);
                }
            }
        }

        // st false cases
        for( auto node : *(splitTreeInfo->vertices) ){
            if( node->criticalType == 3 ){
                int i = node->globalIndex;
                if ( !splitTreeInfo2->tree->hasNode(i) || splitTreeInfo2->tree->getParentGlobalIndex(i) != splitTreeInfo->tree->getParentGlobalIndex(i)
                     || splitTreeInfo2->tree->getNode(i)->criticalType != node->criticalType ){
                    falseNegativesST.push_back(i);
                }
            }

        }

        for( auto node : *(splitTreeInfo2->vertices) ){
            if( node->criticalType == 3 ){
                int i = node->globalIndex;
                if ( !splitTreeInfo->tree->hasNode(i) || splitTreeInfo->tree->getParentGlobalIndex(i) != node->parent->globalIndex ) {
                    falsePositivesST.push_back(i);
                }
            }
        }

        if( falsePositivesJT.size() == 0 && falsePositivesST.size() == 0 && falseNegativesJT.size() == 0 && falseNegativesST.size() == 0 ){
            if( verbose ){
                cout << "it worked :)" << endl;
            }
            keepGoing = false;
        } else {
            if( verbose ){
                cout << "(" << falsePositivesJT.size() << "," << falseNegativesJT.size() << ") (" << falsePositivesST.size() << "," << falseNegativesST.size() << ")" << endl;
            }
            // maps that will store lists of each point needed for segmentation
            unordered_map<int, unique_ptr<vector<int>>> jtBaseSegmentationLists2;
            unordered_map<int, unique_ptr<vector<int>>> stBaseSegmentationLists2;

            // need to repeat this whole business to find the merge tree segmentations.
            vector<int>& jtUSSC2 = *(joinTreeInfo2->USSC);
            vector<int>& jtUSCI2 = *(joinTreeInfo2->USCI);
            vector<int>& stUSSC2 = *(splitTreeInfo2->USSC);
            vector<int>& stUSCI2 = *(splitTreeInfo2->USCI);

            unordered_map<int,int>& jtVTTS2 = *(joinTreeInfo2->visitedTrackerIdToSegmentation);
            unordered_map<int,int>& stVTTS2 = *(splitTreeInfo2->visitedTrackerIdToSegmentation);                        

            // comparisons needed for the following step
            ScalarFieldMore<T> reconIndexIsMore;
            ScalarFieldLess<T> reconIndexIsLess;

            reconIndexIsMore.setSF(cd->reconstructedData);
            reconIndexIsLess.setSF(cd->reconstructedData);

            for( int i = 0; i < numPoints; ++i ){
                vector<int> jtSegmentations;
                vector<int> stSegmentations;

                // bounds from join tree
                int jtSegmentation = joinTreeInfo2->visitedTracker->entries[i];
                if( jtSegmentation == -1 ){

                    // in this case, the join tree segmentation comes from inside of saddles.
                    jtSegmentations.push_back( jtVTTS2[ jtUSCI2[upper_bound( jtUSSC2.begin(), jtUSSC2.end(), i, reconIndexIsLess ) - jtUSSC2.begin() - 1] ] );
                } else if( jtSegmentation == -2 ){

                    for( auto subSegmentation : joinTreeInfo2->visitedTracker->spillOverMarks[ i ] ){
                        jtSegmentations.push_back( jtVTTS2[subSegmentation] );
                    }

                } else { 
                    jtSegmentations.push_back( jtVTTS2[ jtSegmentation ] );
                }

                // bounds from split tree
                int stSegmentation = splitTreeInfo2->visitedTracker->entries[i];
                if( stSegmentation == -1 ){
                    // in this case, the join tree segmentation comes from inside of saddles.
                    stSegmentations.push_back(stVTTS2[ stUSCI2[upper_bound( stUSSC2.begin(), stUSSC2.end(), i, reconIndexIsMore ) - stUSSC2.begin() - 1] ]);
                } else if( stSegmentation == -2 ){
                    for( auto subSegmentation : splitTreeInfo2->visitedTracker->spillOverMarks[ i ] ){
                        stSegmentations.push_back(stVTTS2[subSegmentation]);
                    }
                } else { 
                    stSegmentations.push_back(stVTTS2[stSegmentation]);
                }

                // then add this point to all relevant segmentation lists
                for( auto jtSegmentationElt : jtSegmentations ){
                    if( !jtBaseSegmentationLists2.count(jtSegmentationElt) ){
                        jtBaseSegmentationLists2[jtSegmentationElt] = unique_ptr<vector<int>>(new vector<int>());
                    }
                    jtBaseSegmentationLists2[jtSegmentationElt]->push_back(i);
                }

                for( auto stSegmentationElt : stSegmentations ){
                    if( !stBaseSegmentationLists2.count(stSegmentationElt) ){
                        stBaseSegmentationLists2[stSegmentationElt] = unique_ptr<vector<int>>(new vector<int>());
                    }
                    stBaseSegmentationLists2[stSegmentationElt]->push_back(i);
                }            
            }

            int radius;
            if( step < 3 ){
                radius = 0;
            } else if( step < 6 ){
                radius = 1;
            } else {
                radius = std::min(step - 4,5);
            }

            int numIntervals = lround(pow(2, step+1));

            vector<Region<T>*> tightenRegions;

            for( auto i : falsePositivesJT ){
                Region<T>* tightenRegion = new Region<T>(*(jtBaseSegmentationLists2[i]), cd->reconstructedData);

                if( step < 6 ){
                    Region<T> saddleArea(i, &sf);
                    saddleArea.grow(radius);

                    tightenRegion->mergeInOtherRegion(saddleArea);
                } else {
                    tightenRegion->grow(radius);
                }

                tightenRegions.push_back(tightenRegion);
            }

            for( auto i : falsePositivesST ){
                Region<T>* tightenRegion = new Region<T>(*(stBaseSegmentationLists2[i]), cd->reconstructedData);

                if( step < 6 ){
                    Region<T> saddleArea(i, &sf);
                    saddleArea.grow(radius);

                    tightenRegion->mergeInOtherRegion(saddleArea);
                } else {
                    tightenRegion->grow(radius);
                }

                tightenRegions.push_back(tightenRegion);
            }

            for( auto i : falseNegativesJT ){
                Region<T>* tightenRegion = new Region<T>(*(jtBaseSegmentationLists[i]), &sf);
                tightenRegion->grow(radius);
            }

            for( auto i : falseNegativesST ){
                Region<T>* tightenRegion = new Region<T>(*(stBaseSegmentationLists[i]), &sf);
                tightenRegion->grow(radius);
            }

            for( auto tightenRegion : tightenRegions ){
                sort(tightenRegion->pointList.begin(), tightenRegion->pointList.end(), groundIndexIsLess);

                vector<T> percentileCutoffs(numIntervals+1);

                for( int i = 0; i <= numIntervals; ++i ){
                    percentileCutoffs[i] = sf.getElement( tightenRegion->pointList[ lround( (tightenRegion->pointList.size()-1) * i / numIntervals ) ] );
                }

                for( auto i : tightenRegion->pointList ){

                    int currentInterval = 0;

                    // tighten the error bounds for each point.
                    while( sf.getElement(i) > percentileCutoffs[currentInterval+1] ){
                        ++currentInterval;
                    }

                    bool boundsChanged = false;
                    if( cd->lowerBounds[i] < percentileCutoffs[currentInterval] ){
                        cd->lowerBounds[i] = percentileCutoffs[currentInterval];
                        boundsChanged = true;
                    }

                    if( cd->upperBounds[i] > percentileCutoffs[currentInterval+1] ){
                        cd->upperBounds[i] = percentileCutoffs[currentInterval+1];
                        boundsChanged = true;
                    }

                    if( boundsChanged ){
                        cd->quantizePoint(i);
                    }
                }

                delete tightenRegion;                
            }

            step += 1;

        }

        // if( falsePositives.size() == 0 && falseNegatives.size() == 0 ){
        //     if( verbose ){
        //         cout << "it worked :)" << endl;
        //     }
        //     keepGoing = false;
        // } else {
        //     if( verbose ){
        //         cout << "(" << falsePositives.size() << "," << falseNegatives.size() << ")" << endl;
        //     }

        //     int radius;
        //     if( step < 3 ){
        //         radius = 0;
        //     } else if( step < 6 ){
        //         radius = 1;
        //     } else {
        //         radius = std::min(step - 4,5);
        //     }

        //     int numIntervals = lround(pow(2, step+1));

        //     while( falsePositives.size() > 0 || falseNegatives.size() > 0 ){

        //         int leaf;
        //         bool isFp;

        //         if( falsePositives.size() > 0 ){
        //             leaf = falsePositives.back();
        //             falsePositives.pop_back();
        //             isFp = true;
        //         } else {
        //             leaf = falseNegatives.back();
        //             falseNegatives.pop_back();
        //             isFp = false;
        //         }

        //         Region<T>* tightenRegion = NULL;

        //         if( isFp ){
        //             tightenRegion = new Region<T>((*reconCTISegmentation.CTSegmentationIdToPointList)[(*reconCTISegmentation.leafIdToCTSegmentationId)[leaf]], &sf);

        //             if( step < 6 ){
        //                 Region<T> saddleArea(leaf, &sf);
        //                 saddleArea.grow(radius);

        //                 tightenRegion->mergeInOtherRegion(saddleArea);
        //             } else {
        //                 tightenRegion->grow(radius);
        //             }

        //         } else {
        //             tightenRegion = new Region<T>((*contourTreeInfoSegmentation.CTSegmentationIdToPointList)[(*contourTreeInfoSegmentation.leafIdToCTSegmentationId)[leaf]], &sf);
        //             tightenRegion->grow(radius);
        //         }

        //         // tighten
        //         sort(tightenRegion->pointList.begin(), tightenRegion->pointList.end(), groundIndexIsLess);

        //         vector<T> percentileCutoffs(numIntervals+1);

        //         for( int i = 0; i <= numIntervals; ++i ){
        //             percentileCutoffs[i] = sf.getElement( tightenRegion->pointList[ lround( (tightenRegion->pointList.size()-1) * i / numIntervals ) ] );
        //         }

        //         for( auto i : tightenRegion->pointList ){

        //             int currentInterval = 0;

        //             // tighten the error bounds for each point.
        //             while( sf.getElement(i) > percentileCutoffs[currentInterval+1] ){
        //                 ++currentInterval;
        //             }

        //             bool boundsChanged = false;
        //             if( cd->lowerBounds[i] < percentileCutoffs[currentInterval] ){
        //                 cd->lowerBounds[i] = percentileCutoffs[currentInterval];
        //                 boundsChanged = true;
        //             }

        //             if( cd->upperBounds[i] > percentileCutoffs[currentInterval+1] ){
        //                 cd->upperBounds[i] = percentileCutoffs[currentInterval+1];
        //                 boundsChanged = true;
        //             }

        //             if( boundsChanged ){
        //                 cd->quantizePoint(i);
        //             }
        //         }

        //         delete tightenRegion;
        //         ++step;

        //     }

        //     for( auto fp : falsePositives ){
        //         Region<T>* tightenRegion = NULL;
        //     }
        // }
    }

    auto tightenSplit = std::chrono::high_resolution_clock::now();

    // write to file

    // set up quantization codes and huffman encode

    // set up quantization codes:
    vector<int> quantizationCodes;
    vector<float> losslessStorage32;
    vector<double> losslessStorage64;
    quantizationCodes.reserve(numPoints);

    int quantizationCodeLossless32 = lround(pow( 2, LOSSLESS_PRECISION_32 + 5 )) + 1;
    int quantizationCodeLossless64 = lround(pow( 2, LOSSLESS_PRECISION_64 + 5 )) + 1;

    for( int i = 0; i < numPoints; ++i ){
        if( cd->precisions[i] == LOSSLESS_PRECISION_32 ){
            quantizationCodes.push_back( quantizationCodeLossless32 );
            losslessStorage32.push_back( (float)(sf.getElement(i)) );
        } else if( cd->precisions[i] == LOSSLESS_PRECISION_64 ){
            if( is_same<T,float>() ){
                cd->precisions[i] = LOSSLESS_PRECISION_32;
                quantizationCodes.push_back( quantizationCodeLossless32 );                
                losslessStorage32.push_back( (float)(sf.getElement(i)) );
            } else {
                quantizationCodes.push_back( quantizationCodeLossless64 );                
                losslessStorage64.push_back( sf.getElement(i) );
            }
        } else {
            quantizationCodes.push_back(lround( pow( 2, MAX_PRECISION - cd->precisions[i] ) * cd->quantizations[i] ));
        }
    }

    pair<vector<int>, vector<char>> encodeOutput; 
    encodeOutput = huffmanEncode( quantizationCodes );

    ofstream outf( ( outputFolder + "/values.bytes" ).c_str(), ios::out | ios::binary );
    write_int( outf, sf.size_x );
    write_int( outf, sf.size_y );
    write_int( outf, sf.size_z );
    write_T<T>( outf, xi );
    write_int( outf, losslessStorage32.size() );
    write_int( outf, losslessStorage64.size() );
    write_int( outf, encodeOutput.first.size() ); // will be 0 if we are not using entropy coding. Not a big deal
    write_int( outf, encodeOutput.second.size() ); // same as above.
    write_T<double>( outf, compressorParameter );

    outf.write( (char*)losslessStorage32.data(), losslessStorage32.size()*4 );
    outf.write( (char*)losslessStorage64.data(), losslessStorage64.size()*8 );
    outf.write( (char*)encodeOutput.first.data(), encodeOutput.first.size()*4 );
    outf.write( encodeOutput.second.data(), encodeOutput.second.size() );
    outf.close();

    string suffix;
    
    string cwd = filesystem::current_path();
    int result = chdir(outputFolder.c_str());
    result = system(("tar cvf compressed.tar ./rawData.cmp ./values.bytes" + systemSuffix).c_str());
    result = chdir(cwd.c_str());

    result = system(("rm " + outputFolder + "/compressed.tar.xz" + systemSuffix).c_str());
    result = system(("xz -v9 " + outputFolder + "/compressed.tar" + systemSuffix).c_str());
    suffix = "xz";

    int baseFileSize = filesystem::file_size( filename.c_str() );
    int compressedFileSize = filesystem::file_size( (outputFolder + "/compressed.tar.xz").c_str() );
    double compressionRatio = ((double)baseFileSize) / ((double)compressedFileSize);

    cout << "compression ratio: " << compressionRatio << endl;
    if( outputFilename.length() > 0 ){
        result = system(("cp " + outputFolder + "/compressed.tar.xz " + outputFilename + systemSuffix).c_str());
        cout << "compressed file: " << outputFilename << endl;
    }

    result = system(("rm " + outputFolder + "/intermediateData.dat" + systemSuffix).c_str());
    result = system(("rm " + outputFolder + "/rawData.cmp" + systemSuffix).c_str());
    result = system(("rm " + outputFolder + "/rawData.dat" + systemSuffix).c_str());
    result = system(("rm " + outputFolder + "/values.bytes" + systemSuffix).c_str());

    auto endTime = std::chrono::high_resolution_clock::now();
    double totalTime = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime).count();
    cout << "compression time: " << totalTime << endl;

    results->ratio = compressionRatio;
    results->totalTime = totalTime;

    results->baseTime = std::chrono::duration_cast<std::chrono::duration<double>>(baseCompressionSplit - startTime).count();
    results->gtCtTime = std::chrono::duration_cast<std::chrono::duration<double>>(contourTreeSplit - baseCompressionSplit).count();
    results->errorBoundTime = std::chrono::duration_cast<std::chrono::duration<double>>(errorBoundSplit - contourTreeSplit).count();
    results->growthTime = 0;
    results->averageTightenTime = 0;
    results->numFP = 0;
    results->numFN = 0;
    results->writeToFileTime = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - tightenSplit).count();
    results->numExpanded = 0;

    return results;

}
