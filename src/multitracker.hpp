#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "unionfind.hpp"
#include "unorderedunionfind.hpp"

using namespace std;

class MultiTracker{
    public:
        int * entries;
        unordered_map< int, vector<int> > spillOverMarks;
    public:
        MultiTracker(int size);
        ~MultiTracker();

        // functions with a list based union find.
        void setMark(int position, int mark, UnionFind& uf);
        bool checkMark(int position, int mark, UnionFind& uf );

        // functions with unordered union find.
        void setMark(int position, int mark, UnorderedUnionFind& uf);
        bool checkMark(int position, int mark, UnorderedUnionFind& uf);
        void removeMark(int position, int mark, UnorderedUnionFind& uf);
        void removeMark( int position, int mark );
};