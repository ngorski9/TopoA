#pragma once

#include "evtk.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPointData.h>
using namespace std;

template <typename T> class ScalarField{
    public:
        vector<T> data;
        T dataRange;

        int size_x;
        int size_y;
        int size_z;

        void loadFromVTK(string filename, string arrayName){
            vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            vtkSmartPointer<vtkImageData> imageData = reader->GetOutput();
            vtkSmartPointer<vtkPointData> pointData = imageData->GetPointData();

            vtkSmartPointer<vtkDataArray> dataArray = pointData->GetArray(arrayName.c_str());
            
            if(!dataArray){
                cout << "TopoA: Cannot load array '" + arrayName + "' from VTK file '" + filename + "'" << endl;
                exit(1);
            }

            int dims[3];
            imageData->GetDimensions(dims);
            this->size_x = dims[0];
            this->size_y = dims[1];
            this->size_z = dims[2];
            int numPoints = size_x*size_y*size_z;

            T minValue = INFINITY;
            T maxValue = -INFINITY;

            data.clear();
            data.reserve(numPoints);
            for( int i = 0; i < numPoints; ++i ){
                T val = dataArray->GetTuple1(i);
                data.push_back( val );
                
                if( val < minValue ){
                    minValue = val;
                }

                if( val > maxValue ){
                    maxValue = val;
                }
            }

            dataRange = maxValue - minValue;
        }

        void saveToVTK(string filename){
            saveImageData<T>( filename, size_x, size_y, size_z, data );
        }

        void loadFromDat(string filename, int size_x_, int size_y_, int size_z_){
            size_x = size_x_;
            size_y = size_y_;
            size_z = size_z_;

            int numPoints = size_x * size_y * size_z;
            data = vector<T>(numPoints);

            ifstream in( filename.c_str(), ios::binary );
            in.read( reinterpret_cast<char*>(data.data()), numPoints*sizeof(T) );
            in.close();
        }

        void saveToDat(string filename){
            ofstream out( filename.c_str(), ios::out | ios::binary );
            out.write( (char*)data.data(), data.size()*sizeof(T) );
            out.close();
        }

        void createEmptySF(int size_x_, int size_y_, int size_z_){
            size_x = size_x_;
            size_y = size_y_;
            size_z = size_z_;

            int numPoints = size_x * size_y * size_z;
            data = vector<T>(numPoints);
        }

        ScalarField<float>* downSample(){
            ScalarField<float>* out = new ScalarField<float>();
            out->data = vector<float>( data.begin(), data.end() );
            out->dataRange = dataRange;
            out->size_x = size_x;
            out->size_y = size_y;
            out->size_z = size_z;
            return out;
        }

        ScalarField<double>* upSample(){
            ScalarField<double>* out = new ScalarField<double>();
            out->data = vector<double>( data.begin(), data.end() );
            out->dataRange = dataRange;
            out->size_x = size_x;
            out->size_y = size_y;
            out->size_z = size_z;
            return out;
        }

        T getElement(tuple<int, int, int> coords){
            return data[ coordsToIndex(coords) ];
        }

        T getElement(int x, int y, int z){
            return data[x + y*size_x + z*size_x*size_y];
        }

        T getElement(int idx){
            return data[idx];
        }

        // this breaks the data range feature. Just so you know.
        void setElement(int x, int y, int z, T val){
            data[x + y*size_x + z*size_x*size_y] = val;
        }

        // also breaks the data range feature.
        void setElement(int idx, T val){
            data[idx] = val;
        }

        int coordsToIndex(tuple<int,int,int> coords){
            return get<0>(coords) + get<1>(coords)*size_x + get<2>(coords)*size_x*size_y;
        }

        int coordsToIndex(int x, int y, int z){
            return x + y*this->size_x + z*this->size_x*this->size_y;
        }

        tuple< int, int, int > indexToCoords( int idx ){
            int x = (idx) % size_x;
            int y = ((idx)/size_x)%size_y;
            int z = (idx)/(size_x*size_y);
            return { x, y, z };
        }

        bool less(int x1, int y1, int z1, int x2, int y2, int z2){
            T elt1 = getElement(x1, y1, z1);
            T elt2 = getElement(x2, y2, z2);
            if( elt1 != elt2 ){
                return elt1 < elt2;
            } else {
                return (z1 < z2) || ( z1 == z2 && y1 < y2 ) || ( z1 == z2 && y1 == y2 && x1 < x2 );
            }
        }

        bool less(int idx1, int idx2){
            T elt1 = getElement(idx1);
            T elt2 = getElement(idx2);

            if( elt1 != elt2 ){
                return elt1 < elt2;
            } else {
                return idx1 < idx2;
            }
        }

        bool more(int x1, int y1, int z1, int x2, int y2, int z2){
            T elt1 = getElement(x1, y1, z1);
            T elt2 = getElement(x2, y2, z2);
            if( elt1 != elt2 ){
                return elt1 > elt2;
            } else {
                return (z1 > z2) || ( z1 == z2 && y1 > y2 ) || ( z1 == z2 && y1 == y2 && x1 > x2 );
            }
        }

        int size(){
            return this->size_x * this->size_y * this->size_z;
        }

        T getDataRange(){
            return dataRange;
        }
};