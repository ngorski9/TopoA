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

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

template <typename T> 
Results* compress(string filename, string arrayName, double epsilon_relative, double xi_relative, string outputFilename = "", string outputFolder = ".", string baseCompressor = "SZ3", double compressorParameter = 1, bool verbose = false, int logQuantizeMode = true, int initialPrecision = 0, string baseCompressorFolder = "../../base_compressors", int size_x = -1, int size_y = -1, int size_z = -1){
    Results* results = new Results();

    auto startTime = std::chrono::high_resolution_clock::now();

    ScalarField<T> sf;

    if( size_x == -1 ){
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

    // debug
    // sfFinal.loadFromVTK(filename, arrayName);

    // now for the "progressive" step.

    // set up data structures.
    MultiTracker jtHeapTracker(numPoints);
    MultiTracker stHeapTracker(numPoints);
    MultiTracker jtVisitedTracker(numPoints);
    MultiTracker stVisitedTracker(numPoints);

    UnorderedUnionFind jtGrowthMerge;
    UnorderedUnionFind stGrowthMerge;

    unordered_map< int, ProgressiveCriticalPointWrapper<T>* > jtCPs;
    unordered_map< int, ProgressiveCriticalPointWrapper<T>* > stCPs;

    // to remove things from the heap, we will replace that value with the critical point at the global
    // max or min and pop until it comes out.
    unordered_set<int> jtMinsSet;
    unordered_set<int> stMaxsSet;
    vector<int> jtMinsHeap;
    vector<int> stMaxsHeap;

    jtMinsHeap.reserve( joinTreeInfo->vertices->size() );
    stMaxsHeap.reserve( splitTreeInfo->vertices->size() );

    // comparison functions used for the heaps.
    ScalarFieldLess<T> reconstructedIndexIsLess;
    ScalarFieldMore<T> reconstructedIndexIsMore;

    reconstructedIndexIsLess.setSF(&sfFinal);
    reconstructedIndexIsMore.setSF(&sfFinal);

    // set up the heaps.

    // first, add all of the extrema from the join and split tree (that are critical) to the heaps
    for( auto v : *(joinTreeInfo->vertices) ){

        if( v->criticalType == 0 ){
            jtMinsHeap.push_back( v->globalIndex );
            jtCPs[ v->globalIndex ] = new ProgressiveCriticalPointWrapper<T>(true, v->x, v->y, v->z, &sfFinal);
            jtMinsSet.insert(v->globalIndex);
        }

    }

    for( auto v : *(splitTreeInfo->vertices) ){

        if( v->criticalType == 3 ){
            stMaxsHeap.push_back( v->globalIndex );
            stCPs[ v->globalIndex ] = new ProgressiveCriticalPointWrapper<T>(false, v->x, v->y, v->z, &sfFinal);
            stMaxsSet.insert(v->globalIndex);
        }

    }

    for(int k = 0; k < sf.size_z; ++k){
        for( int j = 0; j < sf.size_y; ++j ){
            for( int i = 0; i < sf.size_x; ++i ){

                tuple<bool, bool> minOrMax = isMinOrMax(&sfFinal, i, j, k);

                if( get<0>(minOrMax) ){
                    int globalIndex = sf.coordsToIndex(i,j,k);
                    if( jtMinsSet.count(globalIndex) == 0 ){
                        jtMinsHeap.push_back( globalIndex );
                        jtCPs[ globalIndex ] = new ProgressiveCriticalPointWrapper<T>(true, i, j, k, &sfFinal);
                        jtMinsSet.insert( globalIndex );
                    }
                } else if( get<1>(minOrMax) ){
                    int globalIndex = sf.coordsToIndex(i,j,k);
                    if( stMaxsSet.count(globalIndex) == 0 ){
                        stMaxsHeap.push_back( globalIndex );
                        stCPs[ globalIndex ] = new ProgressiveCriticalPointWrapper<T>(false, i, j, k, &sfFinal);
                        stMaxsSet.insert( globalIndex );
                    }
                }

            }
        }
    }

    // sort mins from greatest to smallest and maxs from smallest to largest.
    // remember that the heap comparison function works in reverse of normal!
    // (I don't understand why they did it that way either)
    make_heap( jtMinsHeap.begin(), jtMinsHeap.end(), reconstructedIndexIsLess );
    make_heap( stMaxsHeap.begin(), stMaxsHeap.end(), reconstructedIndexIsMore );
            int num = 0;

    int nextPassCode = 0; // used for tracking whether things should be in the heap or not.
    // do the growth step with tightening
    
    // we need to make sure that when we refill the mins heap then everything inside actually ends up getting considered.
    // otherwise it might end up breaking when we get to the bottom.
    int minsTarget = 1;

    int numFP = 0;
    int numFN = 0;
    double FCTime = 0;
    int numExpanded = 0;

    while( jtMinsHeap.size() > minsTarget || stMaxsHeap.size() > 1 ){

        if( verbose ){
            cout << jtMinsHeap.size() << " " << stMaxsHeap.size() << endl;
        }

        // do a little aliasing to simplify the code. Not sure if this will cause a performance hit
        // but I'm guessing it won't.
        bool min = (jtMinsHeap.size() != minsTarget);

        vector<int>& cpHeap = min ? jtMinsHeap : stMaxsHeap;
        unordered_set<int>& heapContents = min ? jtMinsSet : stMaxsSet;
        MultiTracker& heapTracker = min ? jtHeapTracker : stHeapTracker;
        MultiTracker & visitedTracker = min ? jtVisitedTracker : stVisitedTracker;
        UnorderedUnionFind& growthMerge = min ? jtGrowthMerge : stGrowthMerge;
        MergeTreeInfo<T>* mergeTreeInfo = min ? joinTreeInfo : splitTreeInfo;
        unordered_map<int, ProgressiveCriticalPointWrapper<T>*>& cps = min ? jtCPs : stCPs;
        unordered_map<int, unique_ptr<vector<int>>>& baseSegmentationLists = min ? jtBaseSegmentationLists : stBaseSegmentationLists;
        int extremumType = min ? 0 : 3;

        if( !min ){
            minsTarget = 1;
        }

        // standard play for generating the contour tree. Take the node at the front of the heap and grow it.
        int cpIndex = cpHeap[0];

        if( min ){
            pop_heap( cpHeap.begin(), cpHeap.end(), reconstructedIndexIsLess );
        } else {
            pop_heap( cpHeap.begin(), cpHeap.end(), reconstructedIndexIsMore );
        }

        cpHeap.pop_back();

        if( heapContents.count( cpIndex ) ){

            heapContents.erase( heapContents.find( cpIndex ) );
            ++numExpanded;

            growthMerge.add_(cpIndex);
            heapTracker.setMark( cpIndex, cpIndex, growthMerge );
            ProgressiveCriticalPointWrapper<T>* cpw = cps[ cpIndex ];
            shared_ptr<CriticalPoint<T>> cp = cpw->cp;
            tuple<int,int,int> lastTuple = cp->heap->front();

            bool keepGoing = true;

            int falseCase = -1; // -1: nothing is false. 0: false positive. 1: false negative.

            // compute the next edge
            while( keepGoing ){

                if( cp->heap->size <= 0 ){
                    cout << "empty heap" << endl;
                    cout << cpIndex << endl;

                    int numLeft = 0;
                    for( int i = 0; i < numPoints; ++i ){
                        if( heapTracker.checkMark( i, cpIndex, growthMerge ) ){
                            ++numLeft;
                        }
                    }

                    cout << numLeft << endl;

                    exit(1);
                }
                tuple<int,int,int> nextTuple = cp->heap->pop();

                // the visited array is for cleaning up marks, not segmentation. So we add it here.
                cpw->visitedPoints[0]->push_back( sf.coordsToIndex( get<0>(nextTuple), get<1>(nextTuple), get<2>(nextTuple) ) );

                if( (min && sfFinal.less( get<0>(nextTuple), get<1>(nextTuple), get<2>(nextTuple), get<0>(lastTuple), get<1>(lastTuple), get<2>(lastTuple) ))
                ||  (!min && sfFinal.more( get<0>(nextTuple), get<1>(nextTuple), get<2>(nextTuple), get<0>(lastTuple), get<1>(lastTuple), get<2>(lastTuple) ))){
                    // when we hit another saddle, check if we hit a ghost saddle and if it is one then swap out cpw and cp and keep going.
                    // otherwise, terminate and see if we hit a false positive or false negative. (there are no false types for merge trees)
                    int saddleIndex = sfFinal.coordsToIndex(lastTuple);

                    cp->heap->cycledPoints.push_back(saddleIndex);

                    if( cps.find(saddleIndex) != cps.end() ){
                        // if the saddle already exists

                        // find the critical point associated with the position that we just hit.
                        // chain up parents (if applicable). This is only necessary to deal with degenerate saddles.

                        bool degenerate = false;
                        if( cps[saddleIndex]->cp->parent ){
                            degenerate = true;
                        }

                        for( auto c : cps[saddleIndex]->branchChildren ){
                            if( sfFinal.less(cpIndex, c->cp->globalIndex) != min ){
                                degenerate = true;
                            }
                        }

                        ProgressiveCriticalPointWrapper<T>* nextCPW = cps[ saddleIndex ];

                        // if the next point is a ghost point and it is able to expand, then keep
                        // expanding past the next cpw. 
                        // Otherwise, check if it no longer a ghost and check if it is a FN or FP (only possible in degenerate cases).
                        if( degenerate ){
                            // merge the current critical point in with the next one
                            growthMerge.union_( cpIndex, saddleIndex );

                            cp->retire();
                            keepGoing = false;
                            falseCase = 0;
                            cpw->attachBranch(nextCPW);
                        } else if( nextCPW->ghost ){
                            // merge the current critical point in with the next one
                            growthMerge.union_( cpIndex, saddleIndex );

                            if(nextCPW->cp->linkIsVisited(visitedTracker, growthMerge)){
                    
                                cpw->absorbBranch(nextCPW);

                                while( sfFinal.less(saddleIndex, sfFinal.coordsToIndex(cp->heap->front())) != min ){
                                    cp->heap->cycle();
                                }

                                lastTuple = cp->heap->front();

                            } else {
                                keepGoing = false;
                                cp->retire();
                                cpw->attachBranch(nextCPW);
                                if( abs( sfFinal.getElement(cpIndex) - sfFinal.getElement(saddleIndex) ) >= epsilon ){

                                    nextCPW->ghost = false;

                                    if( !( mergeTreeInfo->tree->hasNode(cpIndex) && mergeTreeInfo->tree->getParentGlobalIndex(cpIndex) == saddleIndex ) ){
                                        // false positive
                                        falseCase = 0;
                                    }
                        
                                } else if( mergeTreeInfo->tree->hasNode(cpIndex) && mergeTreeInfo->tree->getParentGlobalIndex(cpIndex) == saddleIndex && mergeTreeInfo->tree->getNode(cpIndex)->criticalType == extremumType){

                                    // false negative
                                    falseCase = 1;
                                }
                            }
                        } else {
                            keepGoing = false;
                            cp->retire();

                            // because this saddle not being a ghost means it was already visited by another extrema closer to it than
                            // the one we just tested that is at least epsilon away from the saddle, we can conclude that the extrema that
                            // we just checked is also at least epsilon away from the saddle.

                            // thus, all that is left is to check if it is a FP or if it is good.
                            if( !( mergeTreeInfo->tree->hasNode(cpIndex) && mergeTreeInfo->tree->getParentGlobalIndex(cpIndex) == saddleIndex && mergeTreeInfo->tree->getNode(cpIndex)->criticalType == extremumType) ){
                                // false positive
                                falseCase = 0;
                            } else {
                                // merge the current critical point in with the next one
                                growthMerge.union_( cpIndex, saddleIndex );                                
                                cpw->attachBranch(nextCPW);                                
                            }
                        }

                    } else {
                        keepGoing = false;
                        cp->retire();

                        // make a new critical point
                        ProgressiveCriticalPointWrapper<T>* nextCPW = new ProgressiveCriticalPointWrapper<T>( min, get<0>(lastTuple), get<1>(lastTuple), get<2>(lastTuple), cpw, &sfFinal );

                        growthMerge.add_(saddleIndex);
                        growthMerge.union_(cpIndex, saddleIndex);
                        cps[ saddleIndex ] = nextCPW;

                        if( abs( sfFinal.getElement(cpIndex) - sfFinal.getElement(saddleIndex) ) < epsilon ){
                            if( mergeTreeInfo->tree->hasNode(cpIndex) && mergeTreeInfo->tree->getParentGlobalIndex(cpIndex) == saddleIndex && mergeTreeInfo->tree->getNode(cpIndex)->criticalType == extremumType ){
                                // false negative
                                falseCase = 1;
                                nextCPW->ghost = true;
                            } else {
                                nextCPW->ghost = true;
                            }
                        } else if( !(mergeTreeInfo->tree->hasNode(cpIndex) && mergeTreeInfo->tree->getParentGlobalIndex(cpIndex) == saddleIndex && mergeTreeInfo->tree->getNode(cpIndex)->criticalType == extremumType) ){
                            // false positive
                            falseCase = 0;
                        }
                    }

                } else {
                    // otherwise, add all neighboring points to the queue
                    // we leave  this for the loop unrolling file.
                    visitedTracker.setMark(sfFinal.coordsToIndex(nextTuple), cpIndex, growthMerge);
                    addNeighborsToHeap<T>(nextTuple, cp->heap, &sfFinal, heapTracker, growthMerge, cpIndex);
                    lastTuple = nextTuple;
                }

            }

            if( verbose ){
                if( falseCase == 0 ){
                    cout << "fp" << endl;
                } else if( falseCase == 1 ){
                    cout << "fn" << endl;
                }
            }

            // in the event of a false case, tighten
            if( falseCase == 0 || falseCase == 1 ){

                if( falseCase == 0 ){
                    ++numFP;
                } else {
                    ++numFN;
                }

                auto FCStart = std::chrono::high_resolution_clock::now();

                Region<T>* tightenRegion = NULL;

                int tighteningStrength = cpw->tighteningStrength;
                // cout << "tightening strength " << tighteningStrength << endl;
                int radius;
                if( tighteningStrength < 3 ){
                    radius = 0;
                } else if( tighteningStrength < 6 ){
                    radius = 1;
                } else {
                    radius = std::min(tighteningStrength - 4,5);
                }

                int numIntervals = lround(pow(2, tighteningStrength+1));
                // first, compute the region
                if( falseCase == 0 ){
                    // false positive
                    tightenRegion = new Region<T>(cpw->visitedPoints, &sf);

                    if( tighteningStrength < 6 ){
                        Region<T> saddleArea(cpIndex, &sf);
                        saddleArea.grow(radius);

                        if( mergeTreeInfo->tree->hasNode(cpIndex) ){

                            Region<T> baseArea(*(baseSegmentationLists[cpIndex]), &sf);
                            baseArea.grow(radius);
                            tightenRegion->mergeInOtherRegion(baseArea);

                        }

                        tightenRegion->mergeInOtherRegion(saddleArea);
                    } else {
                        if( mergeTreeInfo->tree->hasNode(cpIndex) ){
                            Region<T> baseArea(*(baseSegmentationLists[cpIndex]), &sf);
                            tightenRegion->mergeInOtherRegion(baseArea);                            
                        }                        
                        tightenRegion->grow(radius);
                    }

                } else if( falseCase == 1 ){
                    // false negative
                    tightenRegion = new Region<T>(*(baseSegmentationLists[cpIndex]), &sf);
                    tightenRegion->grow(radius);
                }

                // iterate through the points and tighten error bounds according to which interval they land in
                // also detect all overlapping regions so that we can clear them

                vector<int> jtOverlappingRegionsList;
                vector<int> stOverlappingRegionsList;
                unordered_set<int> jtOverlappingRegionsSet;
                unordered_set<int> stOverlappingRegionsSet;

                for( auto i : tightenRegion->pointList ){
                    // find overlapping regions in both the JT and ST
                    int jtHeapMark = jtHeapTracker.entries[i];
                    if( jtHeapMark == -2 ){
                        for( auto mark : jtHeapTracker.spillOverMarks[i] ){
                            if( jtOverlappingRegionsSet.count(mark) == 0 ){
                                jtOverlappingRegionsSet.insert(mark);
                                jtOverlappingRegionsList.push_back(mark);
                            }
                        }
                    } else if( jtHeapMark != -1 && jtOverlappingRegionsSet.count(jtHeapMark) == 0 ){
                        jtOverlappingRegionsSet.insert(jtHeapMark);
                        jtOverlappingRegionsList.push_back(jtHeapMark);
                    }

                    int stHeapMark = stHeapTracker.entries[i];
                    if( stHeapMark == -2 ){
                        for( auto mark : stHeapTracker.spillOverMarks[i] ){
                            if( stOverlappingRegionsSet.count(mark) == 0 ){
                                stOverlappingRegionsSet.insert(mark);
                                stOverlappingRegionsList.push_back(mark);
                            }
                        }
                    } else if( stHeapMark != -1 && stOverlappingRegionsSet.count(stHeapMark) == 0 ){
                        stOverlappingRegionsSet.insert(stHeapMark);
                        stOverlappingRegionsList.push_back(stHeapMark);
                    }
                }


                int currentInterval = 0;

                if( tighteningStrength <= 20 ){

                    // tighten the error bound of points in the region

                    // sort the points in the region in order to derive a percentile based partition
                    sort(tightenRegion->pointList.begin(), tightenRegion->pointList.end(), groundIndexIsLess);

                    vector<T> percentileCutoffs(numIntervals+1);

                    for( int i = 0; i <= numIntervals; ++i ){
                        percentileCutoffs[i] = sf.getElement( tightenRegion->pointList[ lround( (tightenRegion->pointList.size()-1) * i / numIntervals ) ] );
                    }

                    for( auto i : tightenRegion->pointList ){

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
                } else {
                    // otherwise, if we get this far, just store losslessly everything in the region.

                    for( auto i : tightenRegion->pointList ){
                        cd->lowerBounds[i] = 1;
                        cd->upperBounds[i] = 0;
                        cd->quantizePoint(i);
                    }
                }

                ++nextPassCode;

                // grow region by 1 for detecting new extrema
                tightenRegion->grow();

                // go through all regions to discover other regions that they merge into and queue those up as well.
                for( int i = 0; i < jtOverlappingRegionsList.size(); ++i ){
                    int r = jtOverlappingRegionsList[i];

                    // if this branch merges into another branch, then we need to redo that one as well
                    // this needs to continue all the way up.

                    ProgressiveCriticalPointWrapper<T>* branchParent = jtCPs[r]->branchParent;
                    while( branchParent != NULL ){

                        if( branchParent->cp->criticalType == 0 && !jtOverlappingRegionsSet.count(branchParent->cp->globalIndex) ){
                            jtOverlappingRegionsList.push_back(branchParent->cp->globalIndex);
                            jtOverlappingRegionsSet.insert(branchParent->cp->globalIndex);
                        }

                        branchParent = branchParent->branchParent;

                    }

                }

                // go through all regions and unmark all points in those regions
                for( int i = 0; i < jtOverlappingRegionsList.size(); ++i ){
                    int r = jtOverlappingRegionsList[i];

                    // this is going to overvisit points a little bit but that's ok. For now...
                    vector<int>* heapContents = jtCPs[r]->cp->heap->reconstructHeapContents();
                    for( auto pt : *heapContents ){
                        jtHeapTracker.removeMark(pt, r);
                    }
                    delete heapContents;

                    // we only need to redo points visited on the main branch, as
                    // the other branches were not affected, or if they were, this code
                    // will run for them as well.
                    for( auto i : *(jtCPs[r]->visitedPoints[0]) ){
                        jtHeapTracker.removeMark( i, r );
                        jtVisitedTracker.removeMark( i, r );
                    }         

                    for( auto pt : jtCPs[r]->cp->heap->cycledPoints ){
                        jtHeapTracker.removeMark( pt, r );
                    }           
                }

                // for each of the previously visited regions, remove any parent connections,
                // then delete the saddle at the top.
                // if whether the base is an extremum cannot be altered by the tightening
                // (i.e. they do not border the region), reset it and queue it back up.
                for( auto r : jtOverlappingRegionsList ){

                    shared_ptr<CriticalPoint<T>> connectionBase = jtCPs[r]->cp;
                    shared_ptr<CriticalPoint<T>> connectionParent = jtCPs[r]->cp->parent;
                    while( connectionParent ){                        
                        connectionBase->parent = NULL;
                        connectionParent->children.erase( find( connectionParent->children.begin(), connectionParent->children.end(), connectionBase ) );
                        connectionBase = connectionParent;
                        connectionParent = connectionParent->parent;
                    }

                    if( jtCPs[r]->branchParent ){

                        // reset all critical points between the extremum and the saddle
                        jtGrowthMerge.resetBetween_(r, jtCPs[r]->branchParent->cp->globalIndex);

                        ProgressiveCriticalPointWrapper<T>* branchParent = jtCPs[r]->branchParent;

                        if( branchParent->branchChildren.size() == 1 ){
                            int deleteIndex = branchParent->cp->globalIndex;
                            delete jtCPs[deleteIndex];
                            jtCPs.erase( jtCPs.find(deleteIndex) );
                        } else {
                            vector<ProgressiveCriticalPointWrapper<T>*>& parentChildren = branchParent->branchChildren;                            
                            parentChildren.erase( find( parentChildren.begin(), parentChildren.end(), jtCPs[r] ) );
                        }

                        jtCPs[r]->branchParent = NULL;

                    }

                    if( !tightenRegion->contains(r) ){
                        // requeue

                        jtCPs[r]->reset(&sfFinal);
                        jtMinsHeap.push_back(r);
                        jtMinsSet.insert(r);
                        push_heap(jtMinsHeap.begin(), jtMinsHeap.end(), reconstructedIndexIsLess );

                    }
                }

                // now repeat this but for the split tree :(
                // I just copy-pasted it. Can you tell?
                for( int i = 0; i < stOverlappingRegionsList.size(); ++i ){

                    int r = stOverlappingRegionsList[i];   

                    // // if this branch merges into another branch, then we need to redo that one as well
                    // // this needs to continue all the way up.

                    ProgressiveCriticalPointWrapper<T>* branchParent = stCPs[r]->branchParent;
                    while( branchParent != NULL ){

                        if( branchParent->cp->criticalType == 3 && !stOverlappingRegionsSet.count(branchParent->cp->globalIndex) ){
                            stOverlappingRegionsList.push_back(branchParent->cp->globalIndex);
                            stOverlappingRegionsSet.insert(branchParent->cp->globalIndex);
                        }

                        branchParent = branchParent->branchParent;

                    }

                }

                for( int i = 0; i < stOverlappingRegionsList.size(); ++i ){
                    int r = stOverlappingRegionsList[i];   

                    // this is going to overvisit points a little bit but that's ok. For now...
                    vector<int>* heapContents = stCPs[r]->cp->heap->reconstructHeapContents();
                    for( auto pt : *heapContents ){
                        stHeapTracker.removeMark(pt, r);
                    }
                    delete heapContents;     

                    // we only need to redo points visited on the main branch, as
                    // the other branches were not affected, or if they were, this code
                    // will run for them as well.
                    for( auto i : *(stCPs[r]->visitedPoints[0]) ){
                        stHeapTracker.removeMark( i, r );
                        stVisitedTracker.removeMark( i, r );
                    }

                    for( auto pt : stCPs[r]->cp->heap->cycledPoints ){
                        stHeapTracker.removeMark(pt, r);
                    }
                }

                // for each of the previously visited regions, delete the saddle at the top.
                // if whether the base is an extremum cannot be altered by the tightening
                // (i.e. they do not border the region), reset it and queue it back up.
                for( auto r : stOverlappingRegionsList ){

                    shared_ptr<CriticalPoint<T>> connectionBase = stCPs[r]->cp;
                    shared_ptr<CriticalPoint<T>> connectionParent = stCPs[r]->cp->parent;
                    while( connectionParent ){

                        connectionBase->parent = NULL;
                        connectionParent->children.erase( find( connectionParent->children.begin(), connectionParent->children.end(), connectionBase ) );
                        connectionBase = connectionParent;
                        connectionParent = connectionParent->parent;
                    }

                    if( stCPs[r]->branchParent ){

                        stGrowthMerge.resetBetween_( r, stCPs[r]->branchParent->cp->globalIndex );

                        ProgressiveCriticalPointWrapper<T>* branchParent = stCPs[r]->branchParent;

                        if( branchParent->branchChildren.size() == 1 ){
                            int deleteIndex = branchParent->cp->globalIndex;
                            delete stCPs[deleteIndex];
                            stCPs.erase( stCPs.find(deleteIndex) );
                        } else {
                            vector<ProgressiveCriticalPointWrapper<T>*>& parentChildren = branchParent->branchChildren;
                            parentChildren.erase( find( parentChildren.begin(), parentChildren.end(), stCPs[r] ) );
                        }

                        stCPs[r]->branchParent = NULL;

                    }

                    if( !tightenRegion->contains(r) ){
                        // requeue

                        stCPs[r]->reset(&sfFinal);
                        stMaxsHeap.push_back(r);
                        stMaxsSet.insert(r);
                        push_heap(stMaxsHeap.begin(), stMaxsHeap.end(), reconstructedIndexIsMore );

                    }
                }

                // go through the main region and look for any mins or maxes and requeue them.
                int idx = 0;
                for( auto i : tightenRegion->pointList ){
                    ++idx;
                    tuple<int,int,int> coords = sfFinal.indexToCoords(i);
                    tuple<bool, bool> minOrMax = isMinOrMax( &sfFinal, get<0>(coords), get<1>(coords), get<2>(coords) );

                    if( get<0>( minOrMax ) ){

                        jtMinsSet.insert(i);

                        if( jtCPs.find( i ) == jtCPs.end() || jtCPs[i]->cp->criticalType != 0 ){
                            // add new min to the queue
                            jtMinsHeap.push_back(i);
                            push_heap(jtMinsHeap.begin(), jtMinsHeap.end(), reconstructedIndexIsLess);
                            jtCPs[i] = new ProgressiveCriticalPointWrapper<T>( true, get<0>(coords), get<1>(coords), get<2>(coords), &sfFinal );
                        } else {
    
                            // requeue
                            jtCPs[i]->reset(&sfFinal);
                            jtMinsHeap.push_back(i);
                            push_heap(jtMinsHeap.begin(), jtMinsHeap.end(), reconstructedIndexIsLess);
                            if( i == cpIndex ){
                                ++(jtCPs[i]->tighteningStrength);
                            }
                        }

                    } else if( jtCPs.find(i) != jtCPs.end() ){

                        delete jtCPs[i];

                        jtCPs.erase( jtCPs.find(i) );

                        if( jtMinsSet.count(i) ){
                            jtMinsSet.erase( jtMinsSet.find(i) );
                        }

                    }

                    // analogous reverse for split tree
                    if( get<1>( minOrMax ) ){
                        stMaxsSet.insert(i);

                        if( stCPs.find( i ) == stCPs.end() || stCPs[i]->cp->criticalType != 3){
                            // add new min to the queue
                            stMaxsHeap.push_back(i);
                            push_heap(stMaxsHeap.begin(), stMaxsHeap.end(), reconstructedIndexIsMore);
                            stCPs[i] = new ProgressiveCriticalPointWrapper<T>( false, get<0>(coords), get<1>(coords), get<2>(coords), &sfFinal );
                        } else {
                            
                            // requeue
                            stCPs[i]->reset(&sfFinal);
                            stMaxsHeap.push_back(i);
                            push_heap(stMaxsHeap.begin(), stMaxsHeap.end(), reconstructedIndexIsMore);
                            if( i == cpIndex ){
                                ++(stCPs[i]->tighteningStrength);
                            }
                        }

                    } else if( stCPs.find(i) != stCPs.end() ){
                        
                        delete stCPs[i];

                        stCPs.erase( stCPs.find(i) );

                        if( stMaxsSet.count(i) ){
                            stMaxsSet.erase( stMaxsSet.find(i) );
                        }
                    }
                }
                delete tightenRegion;

                auto FCEnd = std::chrono::high_resolution_clock::now();
                FCTime += std::chrono::duration_cast<std::chrono::duration<double>>(FCEnd - FCStart).count();

            } // end of if falseCase == 0 || falseCase == 1 (i.e. there is a false case)
        } // end if cps.find(cpIndex) != cps.end() (i.e. it is a valid node that we are removing from the heap)
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
    results->growthTime = std::chrono::duration_cast<std::chrono::duration<double>>(tightenSplit - errorBoundSplit).count() - FCTime;
    results->averageTightenTime = FCTime / ((double)(numFP + numFN));
    results->numFP = numFP;
    results->numFN = numFN;
    results->writeToFileTime = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - tightenSplit).count();
    results->numExpanded = numExpanded;

    return results;

}
