#include "multitracker.hpp"
#include <iostream>
using namespace std;

MultiTracker::MultiTracker(int size){
    entries = new int[size];
    fill_n(entries, size, -1);
}

MultiTracker::~MultiTracker(){
    delete [] entries;
}

void MultiTracker::setMark(int position, int mark, UnionFind& uf){
    int currentMark = entries[position];
    if( currentMark == -1 ){
        entries[position] = uf.find_(mark);
    }else if( currentMark == -2 ){
        spillOverMarks[position].push_back(uf.find_(mark));
    } else{
        spillOverMarks[position] = vector<int>();
        vector<int>& newVector = spillOverMarks[position];
        
        newVector.reserve(6);
        newVector.push_back(currentMark);
        newVector.push_back(uf.find_(mark));

        entries[position] = -2;
    }
}

bool MultiTracker::checkMark(int position, int mark, UnionFind& uf ){
    int currentMark = entries[position];
    if( currentMark == -1 ){
        return false;
    } else if( currentMark == -2 ){
        for(const int& presentMark : spillOverMarks[position]){
            if( uf.find_(presentMark) == uf.find_(mark) ){
                return true;
            }
        }
        return false;
    } else {
        return uf.find_(mark) == uf.find_(currentMark);
    }
}

void MultiTracker::setMark(int position, int mark, UnorderedUnionFind& uf){
    int currentMark = entries[position];
    if( currentMark == -1 ){
        entries[position] = mark;
    }else if( currentMark == -2 ){
        spillOverMarks[position].push_back(mark);
    } else{
        spillOverMarks[position] = vector<int>();
        vector<int>& newVector = spillOverMarks[position];

        newVector.reserve(6);
        newVector.push_back(currentMark);
        newVector.push_back(mark);

        entries[position] = -2;
    }
}

bool MultiTracker::checkMark(int position, int mark, UnorderedUnionFind& uf ){
    int currentMark = entries[position];
    if( currentMark == -1 ){
        return false;
    } else if( currentMark == -2 ){
        for(const int& presentMark : spillOverMarks[position]){
            if( uf.find_(presentMark) == uf.find_(mark) ){
                return true;
            }
        }
        return false;
    } else {
        return uf.find_(mark) == uf.find_(currentMark);
    }
}

// can lead to undefined behavior!
void MultiTracker::removeMark( int position, int mark, UnorderedUnionFind& uf ){
    int currentMark = entries[position];
    if( currentMark == -2 ){
        for( int i = spillOverMarks[position].size() - 1; i >= 0; --i ){
            if( uf.find_(mark) == uf.find_(spillOverMarks[position][i]) ){
                spillOverMarks[position].erase(spillOverMarks[position].begin() + i);
            }
        }

        if( spillOverMarks[position].size() == 1 ){
            entries[position] = spillOverMarks[position][0];
        } else if( spillOverMarks[position].size() == 0 ){
            entries[position] = -1;
        }
    } else if( currentMark != -1 && uf.find_(mark) == uf.find_(currentMark) ){
        entries[position] = -1;
    }
}

// can lead to undefined behavior!
void MultiTracker::removeMark( int position, int mark ){
    int currentMark = entries[position];
    if( currentMark == -2 ){
        for( int i = spillOverMarks[position].size() - 1; i >= 0; --i ){
            if( mark == spillOverMarks[position][i] ){
                spillOverMarks[position].erase(spillOverMarks[position].begin() + i);
            }
        }

        if( spillOverMarks[position].size() == 1 ){
            entries[position] = spillOverMarks[position][0];
        } else if( spillOverMarks[position].size() == 0 ){
            entries[position] = -1;
        }
    } else if( currentMark != -1 && mark == currentMark ){
        entries[position] = -1;
    }
}