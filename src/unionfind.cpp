#include "unionfind.hpp"
#include <iostream>

UnionFind::UnionFind(int size){
    entries = vector<int>(size);
    for(int i = 0; i < size; ++i){
        entries[i] = i;
    }
}

void UnionFind::union_(int first, int second){
    entries[find_(first)] = second;

}

int UnionFind::find_(int idx){
    int next_idx = entries[idx];
    int temp = -1;

    while( idx != next_idx ){
        temp = next_idx;
        next_idx = entries[next_idx];
        entries[idx] = next_idx;
        idx = temp;
    }

    return idx;
}