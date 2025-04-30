#pragma once

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <unordered_set>
using namespace std;

class UnorderedUnionFind{
    private:
        unordered_map<int,int> entries;
        unordered_map<int,int> originalTreeTracker; // measures exactly which node merges into which node. Used for removing a vertex.
        unordered_map<int,vector<int>> reverseTracking; // measures what merges into each individual node. Used for removing a vertex.
    public:
        void add_(int idx);
        void resetBetween_(int child, int parent);
        void reset_(int idx);
        void union_(int first, int second);
        int find_(int idx);
};