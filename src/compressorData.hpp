#pragma once

#include "scalarfield.hpp"

#include <string>

// the max precision used in a quantization interval
const int MAX_PRECISION = 10;
// max precision + 1 stores the point losslessly in 32 bits. lossless + 2 stores in 64 bits.
const int LOSSLESS_PRECISION_32 = MAX_PRECISION + 1;
const int LOSSLESS_PRECISION_64 = MAX_PRECISION + 2;

template<typename T>
class CompressorData{
public:
    ScalarField<T>* groundData;
    ScalarField<T>* intermediateData;
    ScalarField<T>* reconstructedData;
    vector<T> lowerBounds;
    vector<T> upperBounds;
    vector<int> quantizations;
    vector<int> precisions;
    int numPointsInitialized = 0;
    T xi;

    bool logQuantizeMode;
    int initialPrecision;

    CompressorData( ScalarField<T>* groundData_, ScalarField<T>* intermediateData_, ScalarField<T>* reconstructedData_, T xi_, bool logQuantizeMode_ = true, int initialPrecision_ = 0 ){
                        groundData = groundData_;
                        intermediateData = intermediateData_;
                        reconstructedData = reconstructedData_;
                        xi = xi_;
                        logQuantizeMode = logQuantizeMode_;
                        initialPrecision = initialPrecision_;

                        int numPoints = groundData_->size();
                        lowerBounds.reserve(numPoints);
                        upperBounds.reserve(numPoints);
                        quantizations.reserve(numPoints);
                        precisions.reserve(numPoints);
                    }

    void logQuantize( int idx ){
        // we assume that precision is only set in this function
        // thus, if precision is already equal to max precision, we assume
        // that this function has already taken care of what happens if we hit max precision.
        if( precisions[idx] == LOSSLESS_PRECISION_64 ){
            return;
        }

        T intervalRange = upperBounds[idx] - lowerBounds[idx];

        T lowerBound = lowerBounds[idx] + 0.00001 * intervalRange;
        T upperBound = upperBounds[idx] - 0.00001 * intervalRange;

        int quantizationNumber;
        T predictedValue;

        // quantizationNumber = round( ( groundData->getElement(idx) - intermediateData->getElement(idx) ) * ( pow(2, precisions[idx]-1) ) / xi );

        if( intermediateData->getElement(idx) < lowerBound ){
            quantizationNumber = ceil( ( lowerBound - intermediateData->getElement(idx) ) * ( pow(2, precisions[idx]-1) ) / xi );
        } else if( intermediateData->getElement(idx) > upperBound ){
            quantizationNumber = floor( ( upperBound - intermediateData->getElement(idx) ) * ( pow(2, precisions[idx]-1) ) / xi );
        } else {
            quantizationNumber = 0;
        }

        predictedValue = intermediateData->getElement(idx) + (quantizationNumber) * xi / ( pow(2, precisions[idx]-1) );

        if( groundData->getElement(idx) < lowerBound ){
            lowerBound = groundData->getElement(idx);
        }

        if( groundData->getElement(idx) > upperBound ){
            upperBound = groundData->getElement(idx);
        }

        while( predictedValue < lowerBound || predictedValue > upperBound ){
            ++precisions[idx];

            if( precisions[idx] == LOSSLESS_PRECISION_32 ){
                quantizations[idx] = 0;
                reconstructedData->setElement(idx, (float)(groundData->getElement(idx)));
            } else if( precisions[idx] == LOSSLESS_PRECISION_64){
                reconstructedData->setElement(idx, groundData->getElement(idx));
                return;
            } else {

                // quantizationNumber = round( ( groundData->getElement(idx) - intermediateData->getElement(idx) ) * ( pow(2, precisions[idx]-1) ) / xi );

                if( intermediateData->getElement(idx) < lowerBound ){
                    quantizationNumber = ceil( ( lowerBound - intermediateData->getElement(idx) ) * ( pow(2, precisions[idx]-1) ) / xi );
                } else {
                    quantizationNumber = floor( ( upperBound - intermediateData->getElement(idx) ) * ( pow(2, precisions[idx]-1) ) / xi );
                }

                predictedValue = intermediateData->getElement(idx) + (quantizationNumber) * xi / ( pow(2, precisions[idx]-1) );                
            }

        }

        reconstructedData->setElement(idx, predictedValue);
        quantizations[idx] = quantizationNumber;
    }

    void linearQuantize( int idx ){
        // we assume that precision is only set in this function
        // thus, if precision is already equal to max precision, we assume
        // that this function has already taken care of what happens if we hit max precision.
        if( precisions[idx] == LOSSLESS_PRECISION_64 ){
            return;
        }

        T intervalRange = upperBounds[idx] - lowerBounds[idx];

        T lowerBound = lowerBounds[idx] + 0.00001 * intervalRange;
        T upperBound = upperBounds[idx] - 0.00001 * intervalRange;

        int quantizationNumber;
        T predictedValue;

        // quantizationNumber = round( ( groundData->getElement(idx) - intermediateData->getElement(idx) ) * ( pow(2, precisions[idx]-1) ) / xi );

        if( intermediateData->getElement(idx) < lowerBound ){
            quantizationNumber = ceil( ( lowerBound - intermediateData->getElement(idx) ) * ( pow(2, precisions[idx]-1) ) / xi );
        } else if( intermediateData->getElement(idx) > upperBound ){
            quantizationNumber = floor( ( upperBound - intermediateData->getElement(idx) ) * ( pow(2, precisions[idx]-1) ) / xi );
        } else {
            quantizationNumber = 0;
        }

        predictedValue = intermediateData->getElement(idx) + (quantizationNumber) * xi / ( pow(2, precisions[idx]-1) );

        if( groundData->getElement(idx) < lowerBound ){
            lowerBound = groundData->getElement(idx);
        }

        if( groundData->getElement(idx) > upperBound ){
            upperBound = groundData->getElement(idx);
        }

        while( predictedValue < lowerBound || predictedValue > upperBound ){
            ++precisions[idx];

            if( precisions[idx] == LOSSLESS_PRECISION_32 ){
                quantizations[idx] = 0;
                reconstructedData->setElement(idx, (float)(groundData->getElement(idx)));
            } else if( precisions[idx] == LOSSLESS_PRECISION_64){
                reconstructedData->setElement(idx, groundData->getElement(idx));
                return;
            } else { 
                precisions[idx] = LOSSLESS_PRECISION_32;
            }

        }

        reconstructedData->setElement(idx, predictedValue);
        quantizations[idx] = quantizationNumber;
    }

    void quantizePoint( int idx ){
        if( logQuantizeMode ){
            logQuantize(idx);
        } else {
            linearQuantize(idx);
        }
    }

    void initializeNextPointLossless(){
        reconstructedData->setElement( numPointsInitialized, groundData->getElement( numPointsInitialized ) );
        lowerBounds.push_back(0);
        upperBounds.push_back(0);
        quantizations.push_back(0);
        precisions.push_back(LOSSLESS_PRECISION_64);
        ++numPointsInitialized;
    }

    void initializeNextPointFromBounds(T lowerBound, T upperBound){
        lowerBounds.push_back(lowerBound);
        upperBounds.push_back(upperBound);
        quantizations.push_back(0);
        precisions.push_back(initialPrecision);
        quantizePoint(numPointsInitialized);
        ++numPointsInitialized;
    }

    void initializeNextPointFromBoundsWithoutQuantization(T lowerBound, T upperBound){
        lowerBounds.push_back(lowerBound);
        upperBounds.push_back(upperBound);
        quantizations.push_back(0);
        precisions.push_back(initialPrecision);
        ++numPointsInitialized;        
    }
};