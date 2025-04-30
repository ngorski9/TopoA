#include "evtk.hpp"

using namespace std;

void savePointCloud(string filename, int x[], int y[], int z[], int numPoints){
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

    for(int i = 0; i < numPoints; ++i){
        points->InsertNextPoint( x[i], y[i], z[i] );
        vtkSmartPointer<vtkVertex> nextVertex = vtkSmartPointer<vtkVertex>::New();
        nextVertex->GetPointIds()->SetId(0,i);
        cellArray->InsertNextCell(nextVertex);
    }

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_VERTEX, cellArray);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}

void saveGraph(string filename,int x[], int y[], int z[], int numPoints, vector<pair<int,int>> edges){
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

    for(int i = 0; i < numPoints; ++i){
        points->InsertNextPoint( x[i], y[i], z[i] );
    }

    for( const tuple<int, int>& e : edges ){
        vtkSmartPointer<vtkLine> nextLine  = vtkSmartPointer<vtkLine>::New();
        nextLine->GetPointIds()->SetId(0,get<0>(e));
        nextLine->GetPointIds()->SetId(1,get<1>(e));
        cellArray->InsertNextCell(nextLine);
    }

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_LINE, cellArray);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}
