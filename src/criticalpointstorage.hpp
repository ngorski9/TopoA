#pragma once

#include "scalarfieldheap.hpp"
#include "multitracker.hpp"
#include "unionfind.hpp"

#include <vector>
#include <memory>

template <typename T>
class CriticalPoint{
    public:
        int criticalType;
        bool active; // active means that it is in the list ready for expansion

        int x;
        int y;
        int z;
        int criticalIndex;
        int globalIndex;
        ScalarFieldHeap<T>* heap;
        vector<int> link;

        static int numCP;

        shared_ptr<CriticalPoint<T>> parent = nullptr;
        vector<shared_ptr<CriticalPoint<T>>> children;

        // constructur used for minima
        CriticalPoint(bool min, int criticalIndex_, int x_, int y_, int z_, ScalarField<T>* sf){
            if(min){
                criticalType = 0;
            } else {
                criticalType = 3;
            }
            active = true;
            x = x_;
            y = y_;
            z = z_;
            criticalIndex = criticalIndex_;
            globalIndex = sf->coordsToIndex(x,y,z);
            
            heap = new ScalarFieldHeap<T>(min, sf);
            heap->push( {x,y,z} );
        }

        // constructur used for saddle points
        CriticalPoint(bool min, int criticalIndex_, int x_, int y_, int z_, ScalarFieldHeap<T>* heap_, ScalarField<T>* sf){
            children.reserve(2);
            if(min){
                criticalType = 1;
            } else {
                criticalType = 2;
            }
            active = false;
            x = x_;
            y = y_;
            z = z_;

            criticalIndex = criticalIndex_;
            globalIndex = sf->coordsToIndex(x,y,z);
            heap = heap_;
            findLink(min, sf);
        }

        // constructor used for blank points in the combine step
        CriticalPoint(int criticalType_, int criticalIndex_, int x_, int y_, int z_, ScalarField<T>* sf){
            criticalIndex = criticalIndex_;
            children.reserve(2);
            criticalType = criticalType_;
            x = x_;
            y = y_;
            z = z_;
            globalIndex = sf->coordsToIndex(x,y,z);

            heap = NULL;
        }

        bool linkIsVisited(MultiTracker& visitedTracker, UnionFind& uf){
            for(const int& linkIndex : link){
                if( !visitedTracker.checkMark(linkIndex, criticalIndex, uf) ){
                    return false;
                }
            }

            return true;

        }

        bool linkIsVisited(MultiTracker& visitedTracker, UnorderedUnionFind& uf){
            for(const int& linkIndex : link){
                if( !visitedTracker.checkMark(linkIndex, criticalIndex, uf) ){
                    return false;
                }
            }

            return true;

        }

        void mergeInHeap(shared_ptr<CriticalPoint<T>> other){
            ScalarFieldHeap<T>* mergedHeap = heap->merge(other->heap);

            if( mergedHeap == heap ){
                delete other->heap;
            } else if( mergedHeap == other->heap){
                delete heap;
            }

            heap = mergedHeap;
            other->heap = NULL;
        }

        void mergeInHeapWithoutDeletingOther(shared_ptr<CriticalPoint<T>> other){
            ScalarFieldHeap<T>* mergedHeap = heap->mergeWithoutAltering(other->heap);

            delete heap;
            heap = mergedHeap;

        }

        void retire(){
            heap->retire();
        }

        ~CriticalPoint<T>(){
            delete heap;
        }

    private:

        // called from the constructor.
        void findLink(bool min, ScalarField<T>* sf ){

            globalIndex = sf->coordsToIndex(x,y,z);

            link.reserve(14);

            bool xNotLow = (x != 0);
            bool yNotLow = (y != 0);
            bool zNotLow = (z != 0);

            bool xNotHigh = (x != sf->size_x-1);
            bool yNotHigh = (y != sf->size_y-1);
            bool zNotHigh = (z != sf->size_z-1);

            if( xNotLow && min == sf->less(x-1,y,z,x,y,z) ){
                link.push_back( sf->coordsToIndex(x-1,y,z) );
            }

            if( xNotLow && yNotHigh && min == sf->less(x-1,y+1,z,x,y,z) ){
                link.push_back( sf->coordsToIndex(x-1,y+1,z) );
            }

            if( yNotHigh && min == sf->less(x,y+1,z,x,y,z) ){
                link.push_back( sf->coordsToIndex(x,y+1,z) );
            }

            if( xNotLow && zNotHigh && min == sf->less(x-1,y,z+1,x,y,z) ){
                link.push_back( sf->coordsToIndex(x-1,y,z+1) );
            }

            if( zNotHigh && min == sf->less(x,y,z+1,x,y,z) ){
                link.push_back( sf->coordsToIndex(x,y,z+1) );
            }

            if( yNotHigh && zNotHigh && min == sf->less(x,y+1,z+1,x,y,z) ){
                link.push_back( sf->coordsToIndex(x,y+1,z+1) );
            }

            if( xNotLow && yNotHigh && zNotHigh && min == sf->less(x-1,y+1,z+1,x,y,z) ){
                link.push_back( sf->coordsToIndex(x-1,y+1,z+1) );
            }

            if( xNotHigh && min == sf->less(x+1,y,z,x,y,z) ){
                link.push_back( sf->coordsToIndex(x+1,y,z) );
            }

            if( xNotHigh && yNotLow && min == sf->less(x+1,y-1,z,x,y,z) ){
                link.push_back( sf->coordsToIndex(x+1,y-1,z) );
            }

            if( yNotLow && min == sf->less(x,y-1,z,x,y,z) ){
                link.push_back( sf->coordsToIndex(x,y-1,z) );
            }

            if( xNotHigh && zNotLow && min == sf->less(x+1,y,z-1,x,y,z) ){
                link.push_back( sf->coordsToIndex(x+1,y,z-1) );
            }

            if( zNotLow && min == sf->less(x,y,z-1,x,y,z) ){
                link.push_back( sf->coordsToIndex(x,y,z-1) );
            }

            if( yNotLow && zNotLow && min == sf->less(x,y-1,z-1,x,y,z) ){
                link.push_back( sf->coordsToIndex(x,y-1,z-1) );
            }

            if( xNotHigh && yNotLow && zNotLow && min == sf->less(x+1,y-1,z-1,x,y,z) ){
                link.push_back( sf->coordsToIndex(x+1,y-1,z-1) );
            }

        }


};

template <typename T>
class CriticalPointLess{
    private:
        ScalarField<T>* sf;
    public:
        void setSF(ScalarField<T>* sf_){
            sf = sf_;
        }

        bool operator()(shared_ptr<CriticalPoint<T>>& left, shared_ptr<CriticalPoint<T>>& right){
            T elt1 = sf->getElement(left->globalIndex);
            T elt2 = sf->getElement(right->globalIndex);
            return (elt1 < elt2) || (elt1 == elt2 && left->globalIndex < right->globalIndex);
        }
};

template <typename T>
class CriticalPointMore{
    private:
        ScalarField<T>* sf;
    public:
        void setSF(ScalarField<T>* sf_){
            sf = sf_;
        }

        bool operator()(shared_ptr<CriticalPoint<T>>& left, shared_ptr<CriticalPoint<T>>& right){
            T elt1 = sf->getElement(left->globalIndex);
            T elt2 = sf->getElement(right->globalIndex);
            return (elt1 > elt2) || (elt1 == elt2 && left->globalIndex > right->globalIndex);
        }
};

// represents the branch that is wrapped by the critical point.
template<typename T>
class ProgressiveCriticalPointWrapper{
    public:
        shared_ptr<CriticalPoint<T>> cp = nullptr;
        shared_ptr<CriticalPoint<T>> lastChild = nullptr;
        vector<shared_ptr<vector<int>>> visitedPoints;
        ProgressiveCriticalPointWrapper<T>* branchParent = NULL;
        vector<ProgressiveCriticalPointWrapper<T>*> branchChildren; //branch children.
        int tighteningStrength = 0;
        bool ghost = false;
        int passCode = 0;

        ProgressiveCriticalPointWrapper( shared_ptr<CriticalPoint<T>> cp_ ){
            cp = cp_;
        }

        // used for wrapping extrema
        ProgressiveCriticalPointWrapper( bool min, int x_, int y_, int z_, ScalarField<T>* sf ){
            cp = shared_ptr<CriticalPoint<T>>(new CriticalPoint<T>( min, sf->coordsToIndex(x_, y_, z_), x_, y_, z_, sf ));

            lastChild = cp;
            visitedPoints.push_back( shared_ptr<vector<int>>(new vector<int>()) );
        }

        // used for wrapping saddles
        ProgressiveCriticalPointWrapper(bool min, int x_, int y_, int z_, ProgressiveCriticalPointWrapper<T>* other, ScalarField<T>* sf){
            cp = shared_ptr<CriticalPoint<T>>(new CriticalPoint<T>(min, sf->coordsToIndex(x_, y_, z_), x_, y_, z_, NULL, sf ));

            lastChild = cp;
            other->lastChild->parent = cp;
            cp->children.push_back(other->lastChild);
            other->branchParent = this;

            branchChildren.reserve(2);
            branchChildren.push_back(other);
        }

    private:

        // used for when two branches terminate at the same critical point, but the heaps do not need to combine.
        void attachNodes( ProgressiveCriticalPointWrapper<T>* other ){

            lastChild->parent = other->cp;
            other->cp->children.push_back(lastChild);
            lastChild = other->cp;

            // we leave the ghost node logic for outside of this class
            // because it involves the context of several other nodes.

        }

    public:

        // used for joining heaps and visited points
        // meant to be run after attaching.
        void absorbBranch( ProgressiveCriticalPointWrapper<T>* other ){
            attachNodes(other);

            for( auto c : other->branchChildren ){

                for( auto ptr : c->visitedPoints ){
                    visitedPoints.push_back(ptr);
                }

                cp->heap->mergeMemoryEfficient( c->cp->heap );
            }

            other->branchParent = this;
            branchChildren.push_back(other);

        }

        // used for hooking up two branches that share a common parent, without
        // merging the entire heap and visited points.
        // mostly useful for dealing with degenerate cases where the same saddle can
        // terminate two branches
        void attachBranch( ProgressiveCriticalPointWrapper<T>* other ){
            attachNodes(other);

            branchParent = other;
            other->branchChildren.push_back(this);

            // other->branchParent = this;
            // branchChildren.push_back(other);

        }

        void reset(ScalarField<T>* sf){

            for( auto c : branchChildren ){
                c->branchParent = NULL;
            }
            branchChildren.clear();

            cp->children.clear();
            cp->parent = NULL;
            delete cp->heap;
            cp->heap = new ScalarFieldHeap<T>( (cp->criticalType == 0 || cp->criticalType == 1) , sf);
            cp->heap->push( {cp->x, cp->y, cp->z} );

            branchParent = NULL;
            ghost = false;
            visitedPoints.clear();
            visitedPoints.push_back(shared_ptr<vector<int>>(new vector<int>()));
            lastChild = cp;
        }

        ~ProgressiveCriticalPointWrapper(){
            
            // we don't worry too much about the parent here because if a point is deleted so is its branch parent.
            if( branchParent ){
                branchParent->branchChildren.erase( find(branchParent->branchChildren.begin(), branchParent->branchChildren.end(), this) );
            }
            
            for( auto child : branchChildren ){
                child->branchParent = NULL;
                child->cp->parent = nullptr;
            }
        }
};