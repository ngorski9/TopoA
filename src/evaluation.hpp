#pragma once

#include <cmath>
#include <string>
#include <filesystem>
#include <fstream>
#include <chrono>

#include "criticalpointstorage.hpp"
#include "getCTFull.hpp"
#include "evtk.hpp"
#include "scalarfield.hpp"
#include "results.hpp"

// returns (# FP, #FN, #FC)
template<typename T>
tuple<int, int, int, int> contourTreesMatch( ScalarField<T>* sf1, ScalarField<T>* sf2, double epsilon_, Results* results = NULL ){
    
    ContourTreeInfo cti1 = getContourTree( sf1, epsilon_ );

    auto ct2Start = std::chrono::high_resolution_clock::now();
    ContourTreeInfo cti2 = getContourTree( sf2, epsilon_ );
    auto ct2End = std::chrono::high_resolution_clock::now();

    if( results ){
        results->decompressedCtTime += std::chrono::duration_cast<std::chrono::duration<double>>(ct2End - ct2Start).count();        
        results->unsimplifiedGtCtSize = cti1.unsimplifiedSize;
        results->unsimplifiedDecompressedCtSize = cti2.unsimplifiedSize;
        results->simplifiedCTSize = cti1.contourTreeNodes->size();
    }

    int fp = 0;
    int fn = 0;
    int fc = 0;
    int anyMismatch = 0;

    // only count false cases on the first pass, or else they would get double counted
    for( auto i : *(cti1.contourTreeNodes) ){
        shared_ptr<CriticalPoint<T>> node = cti1.contourTree->getNode(i);
        if( node->criticalType == 0 || node->criticalType == 3 ){
            if ( !cti2.contourTree->hasNode(i) || cti2.contourTree->getParentGlobalIndex(i) != cti1.contourTree->getParentGlobalIndex(i) ){
                ++fn;
                ++anyMismatch;
            } else if ( cti2.contourTree->hasNode(i) && node->criticalType != cti2.contourTree->getNode(i)->criticalType){
                ++fc;
                ++anyMismatch;
            }
        } else if( !cti2.contourTree->hasNode(i) || cti2.contourTree->getParentGlobalIndex(i) != cti1.contourTree->getParentGlobalIndex(i) || node->criticalType != cti2.contourTree->getNode(i)->criticalType ){
            ++anyMismatch;
        }

    }

    for( auto i : *(cti2.contourTreeNodes) ){
        shared_ptr<CriticalPoint<T>> node = cti2.contourTree->getNode(i);
        if( (node->criticalType == 0 || node->criticalType == 3) &&
            ( !cti1.contourTree->hasNode(i) || cti1.contourTree->getParentGlobalIndex(i) != cti2.contourTree->getParentGlobalIndex(i) ) ){
                ++fp;
                ++anyMismatch;
        } else if( !cti1.contourTree->hasNode(i) || cti1.contourTree->getParentGlobalIndex(i) != cti2.contourTree->getParentGlobalIndex(i) ){
            ++anyMismatch;
        }
    }

    return {fp,fn,fc,anyMismatch};
}

template<typename T>
pair<T,double> getMaxErrorAndPSNR( ScalarField<T>* sf1, ScalarField<T>* sf2 ){

    int size = sf1->size();

    T maxError = -INFINITY;
    double mse = 0;

    for( int i = 0; i < size; ++i ){
        T error = fabs( sf1->getElement(i) - sf2->getElement(i) );
        if( error > maxError ){
            maxError = error;
        }
        mse += pow(error, 2);
    }

    mse /= ((double)size);

    double psnr = 10 * ( log10( pow(sf1->getDataRange(), 2) ) - log10( mse ) );

    return {maxError, psnr};

}

template<typename T>
void printEvaluationToConsole( string file1, string array1, string file2, bool secondIsVtk, double xi_, double epsilon_, Results* results, bool checkCorrectness ){

    cout << endl << endl << endl << "evaluating" << endl << endl;

    ScalarField<T>* sf1 = new ScalarField<T>();
    ScalarField<T>* sf2 = new ScalarField<T>();

    sf1->loadFromVTK(file1, array1);

    if( secondIsVtk ){
        sf2->loadFromVTK(file2, "Scalars_");
    } else {
        sf2->loadFromDat(file2, sf1->size_x, sf1->size_y, sf1->size_z);
    }

    if( checkCorrectness ){
        tuple<int,int,int,int> falseCases = contourTreesMatch( sf1, sf2, epsilon_, results );
        cout << "contour tree match: \t";
        if( get<0>(falseCases) == 0 && get<1>(falseCases) == 0 && get<2>(falseCases) == 0 && get<3>(falseCases) == 0){
            cout << "good " << endl;
        } else {
            cout << "bad (" << get<0>(falseCases) << "," << get<1>(falseCases) << "," << get<2>(falseCases) << "," << get<3>(falseCases) << ")" << endl;
        }
    }

    pair<T, double> maxErrorAndPSNR = getMaxErrorAndPSNR(sf1, sf2);
    T maxError = maxErrorAndPSNR.first;
    T xi = sf1->getDataRange()*xi_;
    
    cout << "error bound: \t\t";
    if( maxError <= xi ){
        cout << "good " << endl;
    } else {
        cout << "bad (" << ( maxError / xi ) << ")" << endl;
    }

    cout << endl << "compression ratio: \t" << results->ratio << endl;

    T psnr = maxErrorAndPSNR.second;
    cout << "psnr: \t\t\t" << psnr << endl;

    cout << endl << "compression time: \t" << results->totalTime << endl;
    cout << "decompression time: \t" << results->decompressionTime << endl;

}

template<typename T>
void writeEvaluationToCSV( string csv, string compressor, double compressorParameter, string file1, string array1, string file2, bool secondIsVtk, double xi_, double epsilon_, Results* results, bool iterative, bool iterativeMT, bool logScaleQuantization, int initialPrecision, bool checkCorrectness ){

    bool writeHeader = true;

    if( std::filesystem::exists(csv) ){
        writeHeader = false;
    }

    ScalarField<T>* sf1 = new ScalarField<T>();
    ScalarField<T>* sf2 = new ScalarField<T>();

    sf1->loadFromVTK(file1, array1);

    auto loadStart = std::chrono::high_resolution_clock::now();
    if( secondIsVtk ){
        sf2->loadFromVTK(file2, "Scalars_");
    } else {
        sf2->loadFromDat(file2, sf1->size_x, sf1->size_y, sf1->size_z);
    }
    auto loadEnd = std::chrono::high_resolution_clock::now();

    if( results ){
        results->decompressedCtTime = std::chrono::duration_cast<std::chrono::duration<double>>(loadEnd - loadStart).count();
    }

    tuple<int,int,int,int> falseCases;
    if( checkCorrectness ){
        falseCases = contourTreesMatch( sf1, sf2, epsilon_, results );
    } else {
        falseCases = {-1,-1,-1,-1};
    }

    pair<T, double> maxErrorAndPSNR = getMaxErrorAndPSNR(sf1, sf2);
    T maxError = maxErrorAndPSNR.first;
    T xi = sf1->getDataRange()*xi_;

    string status;
    if( checkCorrectness == (get<0>(falseCases) == 0 && get<1>(falseCases) == 0 && get<2>(falseCases) == 0 && get<3>(falseCases) == 0) && maxError <= xi ){
        status = "good";
    } else {
        status = "bad";
    }

    if( iterative ){
        compressor += " (iterative)";
    }

    if( iterativeMT ){
        compressor += " (iterative MT)";
    }

    if( !logScaleQuantization ){
        compressor += " (linear)";
    }
    if( initialPrecision != 0 ){
        compressor += " (" + to_string(initialPrecision) + ")";
    }

    string output = file1 + "," + compressor + "," + to_string(compressorParameter) + "," + to_string(epsilon_) + "," + to_string(xi_) + "," + status + "," + to_string(results->ratio) + "," 
                          + to_string(maxErrorAndPSNR.second) + "," + "(" + to_string(get<0>(falseCases)) + " " + to_string(get<1>(falseCases)) + " " + to_string(get<2>(falseCases)) + " " + to_string(get<3>(falseCases)) + "),"
                          + to_string(maxError) + "," + to_string(results->totalTime) + "," + to_string(results->decompressionTime) + "," + to_string(results->baseTime) + ","
                          + to_string(results->gtCtTime) + "," + to_string(results->errorBoundTime) + "," + to_string(results->growthTime) + ",(" + to_string(results->numFP) + " "
                          + to_string(results->numFN) + ")," + to_string(results->averageTightenTime) + "," + to_string(results->writeToFileTime) + "," + to_string(results->decompressedCtTime) + ","
                          + to_string(results->unsimplifiedGtCtSize) + "," + to_string(results->unsimplifiedDecompressedCtSize) + "," + to_string(results->simplifiedCTSize) + "," + to_string(results->numExpanded);
    
    cout << output << endl;

    ofstream out;
    out.open(csv, std::ios_base::app);

    if( writeHeader ){
        out << "filename,compressor,parameter,epsilon,xi,status,ratio,PSNR,false cases (final),max error,compression time,decompression time,initial compression time,";
        out << "ground truth contour tree time,error bound time,growth time,false cases (fixed),average fix time,write to file time,decompressed CT time,";
        out << "gt ct size unsimplified,decompressed ct size unsimplified,ct size simplified,vertices expanded" << endl;
    }
    out << output << endl;

}