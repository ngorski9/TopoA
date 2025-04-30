#pragma once
#include <iostream>

void pythonError(int result, bool verbose){
    if( result != 0 ){
        std::cout << "TopoA: Neurcomp failed." << std::endl;
        if( !verbose ){
            std::cout << "\tRerun using -verbose to see the Python stack trace." << std::endl;
        }
        exit(1);
    }
}