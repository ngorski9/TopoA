#pragma once

#include <vector>
#include <unordered_set>
#include "scalarfield.hpp"
using namespace std;

template <typename T>
class ScalarFieldLess{
    private:
        ScalarField<T>* sf;
    public:
        void setSF(ScalarField<T>* sf_){
            sf = sf_;
        }

        bool operator()(const int& left, const int& right){
            T elt1 = sf->getElement(left);
            T elt2 = sf->getElement(right);
            return (elt1 < elt2) || (elt1 == elt2 && left < right);
        }
};

template <typename T>
class ScalarFieldMore{
    private:
        ScalarField<T>* sf;
    public:
        void  setSF(ScalarField<T>* sf_){
            sf = sf_;
        }

        bool operator()(const int& left, const int& right){
            T elt1 = sf->getElement(left);
            T elt2 = sf->getElement(right);
            return (elt1 > elt2) || (elt1 == elt2 && left > right);
        }
};

template <typename T>
class ScalarFieldHeap{
    public:
        bool min = true;

        // these 3 variables will only be used while the heap is currently active
        vector< int >* contents;
        unordered_set<int>* contentsSet;
        unordered_set<int>* inheritedPoints;
        vector<int> cycledPoints;

        vector<int>* uniquePoints = NULL;
        vector<int> removedInheritedPoints;
        vector<ScalarFieldHeap<T>*> heapsMergedIn;

        ScalarField<T>* sf;
        ScalarFieldLess<T> sfl;
        ScalarFieldMore<T> sfm;
    public:
        int size = 0;

        ScalarFieldHeap(bool min_, ScalarField<T>* sf_){
            min  = min_;
            sf = sf_;
            sfl.setSF(sf);
            sfm.setSF(sf);

            contents = new vector<int>();
            contentsSet = new unordered_set<int>();
            inheritedPoints = new unordered_set<int>();
        }

        void cycle(){
            cycledPoints.push_back((*contents)[0]);
            auto pt = pop();
        }

        void push( tuple<int, int, int> elt ){
            int idx = sf->coordsToIndex( get<0>(elt), get<1>(elt), get<2>(elt) );
            if( !contentsSet->count(idx) ){
                contents->push_back(idx);
                contentsSet->insert(idx);

                if(min){
                    push_heap(contents->begin(), contents->end(), sfm);
                } else{
                    push_heap(contents->begin(), contents->end(), sfl);
                }

                size += 1;
            }
        }

        tuple<int, int, int> pop(){
            int index = (*contents)[0];
            if( inheritedPoints->count(index) ){
                removedInheritedPoints.push_back(index);
                inheritedPoints->erase(inheritedPoints->find(index));
            }
            contentsSet->erase(contentsSet->find(index));

            if( min ){
                pop_heap(contents->begin(), contents->end(), sfm);
            } else {
                pop_heap(contents->begin(), contents->end(), sfl);                
            }

            contents->pop_back();

            size -= 1;
            return sf->indexToCoords(index);
        }

        tuple<int, int, int> front(){
            return sf->indexToCoords((*contents)[0]);
        }

        ScalarFieldHeap* merge(ScalarFieldHeap* other){
            ScalarFieldHeap* smaller;
            ScalarFieldHeap* larger;

            if(other->size < size){
                smaller = other;
                larger = this;
            } else {
                smaller = this;
                larger = other;
            }

            larger->contents->reserve( smaller->size + larger->size );
            
            int smallSize = smaller->size;

            for(int i = 0; i < smallSize; ++i){
                if( !larger->contentsSet->count((*(smaller->contents))[i]) ){
                    larger->contents->push_back((*(smaller->contents))[i]);
                    larger->contentsSet->insert((*(smaller->contents))[i]);
                    if( min ){
                        push_heap(larger->contents->begin(),larger->contents->end(), sfm );
                    } else {
                        push_heap(larger->contents->begin(),larger->contents->end(), sfl );
                    }
                    ++larger->size;
                }
            }

            return larger;
        }

        ScalarFieldHeap* clone(){
            ScalarFieldHeap<T>* returnHeap = new ScalarFieldHeap<T>(min, sf);

            returnHeap->contents->reserve(contents->size());
            
            for( auto c : *contents ){
                returnHeap->contents->push_back(c);
                returnHeap->contentsSet->insert(c);
            }

            returnHeap->size = size;

            return returnHeap;
        }

        vector<int>* reconstructHeapContents(){
            unordered_set<int> doNotAdd;
            unordered_set<int>* heapContentsSet = new unordered_set<int>();
            vector<int>* returnVector = new vector<int>();
            int size = 0; // just because we have to
            mergeMemoryEfficientHelper(returnVector, heapContentsSet, doNotAdd, size);
            delete heapContentsSet;

            return returnVector;
        }

        void mergeMemoryEfficient(ScalarFieldHeap* other){
            unordered_set<int> doNotAdd;
            other->mergeMemoryEfficientHelper(contents, contentsSet, doNotAdd, size, inheritedPoints );
            heapsMergedIn.push_back(other);
        }

        void mergeMemoryEfficientHelper(vector<int>* heapContentsList, unordered_set<int>* heapContentsSet, unordered_set<int>& doNotAdd, int& size, unordered_set<int>* inheritedPointsSet = NULL ){

            for( auto i : *uniquePoints ){
                if( !heapContentsSet->count(i) && !doNotAdd.count(i) ){
                    heapContentsSet->insert(i);
                    heapContentsList->push_back(i);

                    if( inheritedPointsSet ){
                        inheritedPointsSet->insert(i);
                    }

                    if( min ){
                        push_heap(heapContentsList->begin(),heapContentsList->end(), sfm );
                    } else {
                        push_heap(heapContentsList->begin(),heapContentsList->end(), sfl );
                    }
                    ++size;
                }
            }

            vector<int> addedToDoNotAdd;
            for( auto i : removedInheritedPoints ){
                if( !doNotAdd.count(i) ){
                    addedToDoNotAdd.push_back(i);
                    doNotAdd.insert(i);
                }
            }

            for( auto h : heapsMergedIn ){
                h->mergeMemoryEfficientHelper(heapContentsList, heapContentsSet, doNotAdd, size, inheritedPointsSet);
            }

            for( auto i : addedToDoNotAdd ){
                doNotAdd.erase( doNotAdd.find(i) );
            }
        }

        ScalarFieldHeap* mergeWithoutAltering(ScalarFieldHeap* other){
            ScalarFieldHeap* smaller;
            ScalarFieldHeap* larger;

            if(other->size < size){
                smaller = other;
                larger = this;
            } else {
                smaller = this;
                larger = other;
            }

            ScalarFieldHeap* returnHeap = larger->clone();

            returnHeap->contents->reserve( smaller->size + larger->size );
            
            int smallSize = smaller->size;

            for(int i = 0; i < smallSize; ++i){
                if( !returnHeap->contentsSet->count((*(smaller->contents))[i]) ){
                    returnHeap->contents->push_back((*(smaller->contents))[i]);
                    returnHeap->contentsSet->insert((*(smaller->contents))[i]);
                    
                    if( min ){
                        push_heap(returnHeap->contents->begin(),returnHeap->contents->end(), sfm );
                    } else {
                        push_heap(returnHeap->contents->begin(),returnHeap->contents->end(), sfl );
                    }

                    ++returnHeap->size;
                }
            }

            return returnHeap;
        }

        void retire(){
            uniquePoints = new vector<int>();
            for( auto c : *contents ){
                if( !inheritedPoints->count(c) ){
                    uniquePoints->push_back(c);
                }
            }

            delete contents;
            delete contentsSet;
            delete inheritedPoints;

            contents = NULL;
            contentsSet = NULL;
            inheritedPoints = NULL;
        }

        ~ScalarFieldHeap(){
            delete contents;
            delete contentsSet;
            delete inheritedPoints;
            delete uniquePoints;
        }

};