#pragma once

#include <fstream>

void write_int(ofstream& out, int i){
    out.write( (char*)(&i), sizeof(int) );
}

template<typename T>
void write_T(ofstream& out, T t){
    out.write( (char*)(&t), sizeof(T) );
}