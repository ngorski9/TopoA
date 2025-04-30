#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertex.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

using namespace std;

void savePointCloud(string filename, int x[], int y[], int z[], int numPoints);

void saveGraph(string filename,int x[], int y[], int z[], int numPoints, vector<pair<int,int>> edges);

template<typename T, typename U>
void saveImageDataHelper( string filename, int size_x, int size_y, int size_z, vector<T> data ){

    int numPoints = size_x * size_y * size_z;

    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();

    image->SetDimensions(size_x, size_y, size_z);

    vtkSmartPointer<U> arr = vtkSmartPointer<U>::New();

    arr->SetNumberOfComponents(1);
    arr->SetNumberOfTuples( numPoints );

    for (int i = 0; i < numPoints; ++i) {
        arr->SetTuple1(i, data[i]);
    }

    image->GetPointData()->AddArray(arr);
    arr->SetName("Scalars_");


    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer->SetFileName(filename.c_str());
    writer->SetInputData(image);
    writer->Write();
}

template<typename T>
void saveImageData( string filename, int size_x, int size_y, int size_z, vector<T> data ){
    if( is_same< T, double >::value ){
        saveImageDataHelper<T, vtkDoubleArray>( filename, size_x, size_y, size_z, data );
    } else if( is_same< T, float >::value ){
        saveImageDataHelper<T, vtkFloatArray>( filename, size_x, size_y, size_z, data );
    } else if( is_same< T, int >::value ){
        saveImageDataHelper<T, vtkIntArray>( filename, size_x, size_y, size_z, data );
    } else{
        cout << "you are trying to export to a type of data that is not yet supported" << endl;
        exit(1);
    }
}