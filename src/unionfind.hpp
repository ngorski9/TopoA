#pragma once

#include <vector>
using namespace std;

class UnionFind{
    private:

    public:
        vector<int> entries;    
        UnionFind(int size);

        void union_(int first, int second);
        int find_(int idx);

};