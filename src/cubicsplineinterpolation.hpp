#pragma once

#include <string>
#include <algorithm>

#include "cubicsplinecoefficients.hpp"
#include "scalarfield.hpp"

using namespace std;

int getCompressedDimensionSize( int base_size, int factor ){
    if( base_size % factor == 0 ){
        return base_size / factor + 1;
    } else {
        return base_size / factor + 2;
    }
}

template <typename T>
void csiCompressScalarField(ScalarField<T>& sf, string compressedName, int factor){

    ScalarField<T> outputField;

    int x_size_compressed = getCompressedDimensionSize(sf.size_x, factor);
    int y_size_compressed = getCompressedDimensionSize(sf.size_y, factor);
    int z_size_compressed = getCompressedDimensionSize(sf.size_z, factor);

    outputField.createEmptySF( x_size_compressed, y_size_compressed, z_size_compressed );
    for( int k = 0; k < outputField.size_z; ++k ){
        for( int j = 0; j < outputField.size_y; ++j ){
            for( int i = 0; i < outputField.size_x; ++i ){

                outputField.setElement( i, j, k, sf.getElement(min(i*factor, sf.size_x - 1), min(j*factor, sf.size_y - 1), min(k*factor, sf.size_z - 1) ) );

            }
        }
    }

    outputField.saveToDat(compressedName);
}

template <typename T>
void csiCompressVTK(string filename, string arrayName, string compressedName, int factor){
    ScalarField<T> sf;
    sf.loadFromVTK(filename, arrayName);
    csiCompressScalarField<T>(sf, compressedName, factor);
}

template <typename T>
void csiCompressDat(string filename, int x, int y, int z, string compressedName, int factor){
    ScalarField<T> sf;
    sf.loadFromDat(filename, x, y, z);
    csiCompressScalarField<T>(sf, compressedName, factor);
}

template<typename T>
ScalarField<T>* csiDecompress(string filename, int size_x, int size_y, int size_z, int factor){

    ScalarField<T> compressed;
    compressed.loadFromDat(filename, getCompressedDimensionSize(size_x, factor), getCompressedDimensionSize(size_y, factor), getCompressedDimensionSize(size_z, factor));

    ScalarField<T>* sf = new ScalarField<T>();
    sf->createEmptySF(size_x, size_y, size_z);

    tuple<vector<double>,vector<double>,vector<double>,vector<double>> coefs = splineCoefficients(factor);

    // interpolate in the x direction

    int x_tail_size = (size_x % factor) - 1;
    if( x_tail_size <= 0 ){
        x_tail_size += factor;
    }

    tuple<vector<double>,vector<double>,vector<double>,vector<double>> x_tail_coefs = splineCoefficients(factor, x_tail_size);

    for( int k = 0; k < compressed.size_z; ++k ){
        for( int j = 0; j < compressed.size_y; ++j ){

            int k_ = min(k*factor,sf->size_z-1);
            int j_ = min(j*factor,sf->size_y-1);

            int i0 = sf->coordsToIndex(0, j_, k_);
            double y0 = compressed.getElement(0,j,k);
            double y1 = compressed.getElement(1,j,k);
            double y2 = compressed.getElement(2,j,k);
            double y3 = compressed.getElement(3,j,k);

            // interpolate first two intervals
            for( int i = 0; i < 2*factor; ++i ){
                double elt = get<0>(coefs)[i] * y0 + get<1>(coefs)[i] * y1 + get<2>(coefs)[i] * y2 + get<3>(coefs)[i] * y3;
                sf->setElement( i0 + i, elt );
            }

            // interpolate all intervals except first 2 and last 2
            for( int interval = 2; interval < compressed.size_x - 3; ++interval ){

                i0 = sf->coordsToIndex(interval*factor, j_, k_);
                y0 = y1;
                y1 = y2;
                y2 = y3;
                y3 = compressed.getElement(interval+2,j,k);

                for( int i = 0; i < factor; ++i ){
                    sf->setElement( i0+i, get<0>(coefs)[i + factor] * y0 + get<1>(coefs)[i + factor] * y1 + get<2>(coefs)[i + factor] * y2 + get<3>(coefs)[i + factor] * y3);
                }

            }

            // handle the last two intervals
            i0 = sf->coordsToIndex( sf->size_x - factor - x_tail_size - 1 , j_, k_ );
            y0 = y1;
            y1 = y2;
            y2 = y3;
            y3 = compressed.getElement( compressed.size_x-1, j, k );

            for( int i = 0; i < factor + x_tail_size + 1; ++i ){
                sf->setElement( i0+i, get<0>(x_tail_coefs)[i+factor] * y0 + get<1>(x_tail_coefs)[i+factor] * y1 + get<2>(x_tail_coefs)[i + factor] * y2 + get<3>(x_tail_coefs)[i + factor] * y3 );
            }
        }
    }

    // interpolate in the y direction

    int y_tail_size = (size_y % factor) - 1;
    if( y_tail_size <= 0 ){
        y_tail_size += factor;
    }

    tuple<vector<double>,vector<double>,vector<double>,vector<double>> y_tail_coefs = splineCoefficients(factor, y_tail_size);
    int yShift = sf->size_x;

    for( int k = 0; k < compressed.size_z; ++k ){
        for( int i = 0; i < sf->size_x; ++i ){
            int k_ = min(k*factor,sf->size_z-1);

            int j0 = sf->coordsToIndex(i, 0, k_);
            T y0 = sf->getElement(i,0,k_);
            T y1 = sf->getElement(i,factor,k_);
            T y2 = sf->getElement(i,2*factor,k_);
            T y3 = sf->getElement(i,3*factor,k_);

            // interpolate first two intervals
            for( int j = 0; j < 2*factor; ++j ){
                sf->setElement( j0 + j*yShift, get<0>(coefs)[j] * y0 + get<1>(coefs)[j] * y1 + get<2>(coefs)[j] * y2 + get<3>(coefs)[j] * y3);
            }

            // interpolate all intervals except first 2 and last 2
            for( int interval = 2; interval < compressed.size_y - 3; ++interval ){

                j0 = sf->coordsToIndex(i, interval*factor, k_);
                y0 = y1;
                y1 = y2;
                y2 = y3;
                y3 = sf->getElement(i, factor*(interval+2),k_);

                for( int j = 0; j < factor; ++j ){
                    sf->setElement( j0 + j*yShift, get<0>(coefs)[j + factor] * y0 + get<1>(coefs)[j + factor] * y1 + get<2>(coefs)[j + factor] * y2 + get<3>(coefs)[j + factor] * y3);
                }

            }

            // handle the last two intervals
            j0 = sf->coordsToIndex( i, sf->size_y - factor - y_tail_size - 1 , k_ );
            y0 = y1;
            y1 = y2;
            y2 = y3;
            y3 = sf->getElement( i, sf->size_y-1, k_ );

            for( int j = 0; j < factor + y_tail_size + 1; ++j ){
                sf->setElement( j0 + j*yShift, get<0>(y_tail_coefs)[j+factor] * y0 + get<1>(y_tail_coefs)[j+factor] * y1 + get<2>(y_tail_coefs)[j + factor] * y2 + get<3>(y_tail_coefs)[j + factor] * y3 );
            }
        }
    }

    // interpolate in the z direction

    int z_tail_size = (size_z % factor) - 1;
    if( z_tail_size <= 0 ){
        z_tail_size += factor;
    }

    tuple<vector<double>,vector<double>,vector<double>,vector<double>> z_tail_coefs = splineCoefficients(factor, z_tail_size);
    int zShift = sf->size_x * sf->size_y;

    for( int j = 0; j < sf->size_y; ++j ){
        for( int i = 0; i < sf->size_x; ++i ){

            int k0 = sf->coordsToIndex(i, j, 0);
            T y0 = sf->getElement(i,j,0);
            T y1 = sf->getElement(i,j,factor);
            T y2 = sf->getElement(i,j,2*factor);
            T y3 = sf->getElement(i,j,3*factor);

            // interpolate first two intervals
            for( int k = 0; k < 2*factor; ++k ){
                sf->setElement( k0 + k*zShift, get<0>(coefs)[k] * y0 + get<1>(coefs)[k] * y1 + get<2>(coefs)[k] * y2 + get<3>(coefs)[k] * y3);
            }

            // interpolate all intervals except first 2 and last 2
            for( int interval = 2; interval < compressed.size_z - 3; ++interval ){

                k0 = sf->coordsToIndex(i, j, interval*factor);

                y0 = y1;
                y1 = y2;
                y2 = y3;

                y3 = sf->getElement(i, j, factor*(interval+2) );
                for( int k = 0; k < factor; ++k ){
                    sf->setElement( k0 + k*zShift, get<0>(coefs)[k + factor] * y0 + get<1>(coefs)[k + factor] * y1 + get<2>(coefs)[k + factor] * y2 + get<3>(coefs)[k + factor] * y3);
                }

            }

            // handle the last two intervals
            k0 = sf->coordsToIndex( i, j, sf->size_z - factor - z_tail_size - 1 );
            y0 = y1;
            y1 = y2;
            y2 = y3;
            y3 = sf->getElement( i, j, sf->size_z-1 );

            for( int k = 0; k < factor + z_tail_size + 1; ++k ){
                sf->setElement( k0 + k*zShift, get<0>(z_tail_coefs)[k+factor] * y0 + get<1>(z_tail_coefs)[k+factor] * y1 + get<2>(z_tail_coefs)[k + factor] * y2 + get<3>(z_tail_coefs)[k + factor] * y3 );
            }

        }
    }

    return sf;

}