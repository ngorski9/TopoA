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
#include "getCTFullWithSegmentation.hpp"
#include "cubicsplineinterpolation.hpp"
#include "results.hpp"
#include "write_utils.hpp"
#include "python_error.hpp"

using namespace std;

template <typename T> 
Results* compressIterative(string filename, string arrayName, double epsilon_relative, double xi_relative, string outputFilename = "compressed", string outputFolder = ".", string baseCompressor = "SZ3", double compressorParameter = 1, bool verbose = false, bool logQuantizeMode = true, int initialPrecision = 0, string baseCompressorFolder = "../../base_compressors", int size_x = -1, int size_y = -1, int size_z = -1){
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

    ContourTreeInfoWithSegmentation<T> contourTreeInfoSegmentation = getContourTreeWithSegmentation(&sf, epsilon_relative, true, xi, verbose);
    ContourTreeInfo<T>& contourTreeInfo = *(contourTreeInfoSegmentation.contourTreeInfo);

    MergeTreeInfo<T>* joinTreeInfo = contourTreeInfo.joinTreeInfo;
    MergeTreeInfo<T>* splitTreeInfo = contourTreeInfo.splitTreeInfo;

    auto contourTreeSplit = std::chrono::high_resolution_clock::now();

    // create empty scalar field that will store the final values.
    ScalarField<T> sfFinal = ScalarField<T>();
    sfFinal.createEmptySF(sf.size_x, sf.size_y, sf.size_z);

    // upper bound, lower bound, quantization, precision storage are abstracted
    unique_ptr<CompressorData<T>> cd(new CompressorData<T>(&sf, &sfIntermediate, &sfFinal, xi, logQuantizeMode, initialPrecision));

    // comparisons needed for the following step
    ScalarFieldMore<T> groundIndexIsMore;
    ScalarFieldLess<T> groundIndexIsLess;

    groundIndexIsMore.setSF(&sf);
    groundIndexIsLess.setSF(&sf);

    for( int i = 0; i < numPoints; ++i ){

        if( joinTreeInfo->tree->hasNode(i) || splitTreeInfo->tree->hasNode(i) ){
            cd->initializeNextPointLossless();
        } else {
            T lowerBound = (*contourTreeInfoSegmentation.bounds)[i].first;
            T upperBound = (*contourTreeInfoSegmentation.bounds)[i].second;

            cd->initializeNextPointFromBounds(lowerBound, upperBound);            
        }

    }

    if( verbose ){
        cout << "set initial error bounds" << endl;
    }

    auto errorBoundSplit = std::chrono::high_resolution_clock::now();

    // iterative tightening

    bool keepGoing = true;
    int step = 0;
    while( keepGoing ){
        ContourTreeInfoWithSegmentation<T> reconCTISegmentation = getContourTreeWithSegmentation(cd->reconstructedData, epsilon_relative, false, xi, verbose);
        ContourTreeInfo<T>& contourTreeInfo2 = *(reconCTISegmentation.contourTreeInfo);

        // check for fp and fn
        vector<int> falsePositives;
        vector<int> falseNegatives;

        // only count false types on the first pass, or else they would get double counted
        for( auto i : *(contourTreeInfo.contourTreeNodes) ){
            shared_ptr<CriticalPoint<T>> node = contourTreeInfo.contourTree->getNode(i);
            if( node->criticalType == 0 || node->criticalType == 3 ){
                if ( !contourTreeInfo2.contourTree->hasNode(i) || contourTreeInfo2.contourTree->getParentGlobalIndex(i) != contourTreeInfo.contourTree->getParentGlobalIndex(i)
                     || contourTreeInfo2.contourTree->getNode(i)->criticalType != node->criticalType ){
                    falseNegatives.push_back(node->globalIndex);
                }
            }

        }

        for( auto i : *(contourTreeInfo2.contourTreeNodes) ){
            shared_ptr<CriticalPoint<T>> node = contourTreeInfo2.contourTree->getNode(i);
            if( node->criticalType == 0 || node->criticalType == 3 ){
                if ( !contourTreeInfo.contourTree->hasNode(i) || contourTreeInfo.contourTree->getParentGlobalIndex(i) != contourTreeInfo2.contourTree->getParentGlobalIndex(i) ) {
                    falsePositives.push_back(node->globalIndex);
                }
            }
        }

        if( falsePositives.size() == 0 && falseNegatives.size() == 0 ){
            if( verbose ){
                cout << "it worked :)" << endl;
            }
            keepGoing = false;
        } else {
            if( verbose ){
                cout << "(" << falsePositives.size() << "," << falseNegatives.size() << ")" << endl;
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

            while( falsePositives.size() > 0 || falseNegatives.size() > 0 ){

                int leaf;
                bool isFp;

                if( falsePositives.size() > 0 ){
                    leaf = falsePositives.back();
                    falsePositives.pop_back();
                    isFp = true;
                } else {
                    leaf = falseNegatives.back();
                    falseNegatives.pop_back();
                    isFp = false;
                }

                Region<T>* tightenRegion = NULL;

                if( isFp ){
                    tightenRegion = new Region<T>((*reconCTISegmentation.CTSegmentationIdToPointList)[(*reconCTISegmentation.leafIdToCTSegmentationId)[leaf]], &sf);

                    if( step < 6 ){
                        Region<T> saddleArea(leaf, &sf);
                        saddleArea.grow(radius);

                        tightenRegion->mergeInOtherRegion(saddleArea);
                    } else {
                        tightenRegion->grow(radius);
                    }

                } else {
                    tightenRegion = new Region<T>((*contourTreeInfoSegmentation.CTSegmentationIdToPointList)[(*contourTreeInfoSegmentation.leafIdToCTSegmentationId)[leaf]], &sf);
                    tightenRegion->grow(radius);
                }

                // tighten
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
                ++step;

            }
        }
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
    int compressedFileSize = filesystem::file_size( (outputFolder + "/compressed.tar." + suffix).c_str() );
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
