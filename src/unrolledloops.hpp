#pragma once

#include <tuple>
#include <iostream>
#include "scalarfield.hpp"
#include "multitracker.hpp"
#include "unionfind.hpp"
#include "unorderedunionfind.hpp"

using namespace std;

template<typename T>
inline void addNeighborsToHeap(tuple<int, int, int>& point, ScalarFieldHeap<T>* heap, ScalarField<T>* sf, MultiTracker& visitedTracker, UnionFind& components, int component){
    int x = get<0>(point);
    int y = get<1>(point);
    int z = get<2>(point);

    bool xNotLow = (x != 0);
    bool yNotLow = (y != 0);
    bool zNotLow = (z != 0);

    bool xNotHigh = (x != sf->size_x-1);
    bool yNotHigh = (y != sf->size_y-1);
    bool zNotHigh = (z != sf->size_z-1);

    if( xNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x-1,y,z), component, components )){
        heap->push({x-1,y,z});
        visitedTracker.setMark( sf->coordsToIndex(x-1,y,z), component, components );
    }

    if( xNotLow && yNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x-1,y+1,z), component, components ) ){
        heap->push({x-1,y+1,z});
        visitedTracker.setMark( sf->coordsToIndex(x-1,y+1,z), component, components );            
    }

    if( yNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x,y+1,z), component, components )){
        heap->push({x,y+1,z});
        visitedTracker.setMark( sf->coordsToIndex(x,y+1,z), component, components );
    }

    if( xNotLow && zNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x-1,y,z+1), component, components )){
        heap->push({x-1,y,z+1});
        visitedTracker.setMark( sf->coordsToIndex(x-1,y,z+1), component, components );
    }

    if( zNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x,y,z+1), component, components )){
        heap->push({x,y,z+1});
        visitedTracker.setMark( sf->coordsToIndex(x,y,z+1), component, components );
    }

    if( yNotHigh && zNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x,y+1,z+1), component, components )){
        heap->push({x,y+1,z+1});
        visitedTracker.setMark( sf->coordsToIndex(x,y+1,z+1), component, components );
    }

    if( xNotLow && yNotHigh && zNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x-1,y+1,z+1), component, components )){
        heap->push({x-1,y+1,z+1});
        visitedTracker.setMark( sf->coordsToIndex(x-1,y+1,z+1), component, components );
    }

    if( xNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x+1,y,z), component, components )){
        heap->push({x+1,y,z});
        visitedTracker.setMark( sf->coordsToIndex(x+1,y,z), component, components );
    }

    if( xNotHigh && yNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x+1,y-1,z), component, components )){
        heap->push({x+1,y-1,z});
        visitedTracker.setMark( sf->coordsToIndex(x+1,y-1,z), component, components );
    }

    if( yNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x,y-1,z), component, components )){
        heap->push({x,y-1,z});
        visitedTracker.setMark( sf->coordsToIndex(x,y-1,z), component, components );
    }

    if( xNotHigh && zNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x+1,y,z-1), component, components )){
        heap->push({x+1,y,z-1});
        visitedTracker.setMark( sf->coordsToIndex(x+1,y,z-1), component, components );
    }

    if( zNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x,y,z-1), component, components )){
        heap->push({x,y,z-1});
        visitedTracker.setMark( sf->coordsToIndex(x,y,z-1), component, components );
    }

    if( yNotLow && zNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x,y-1,z-1), component, components )){
        heap->push({x,y-1,z-1});
        visitedTracker.setMark( sf->coordsToIndex(x,y-1,z-1), component, components );
    }

    if( xNotHigh && yNotLow && zNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x+1,y-1,z-1), component, components )){
        heap->push({x+1,y-1,z-1});
        visitedTracker.setMark( sf->coordsToIndex(x+1,y-1,z-1), component, components );
    }
}

template<typename T>
inline void addNeighborsToHeap(tuple<int, int, int>& point, ScalarFieldHeap<T>* heap, ScalarField<T>* sf, MultiTracker& visitedTracker, UnorderedUnionFind& components, int component){
    int x = get<0>(point);
    int y = get<1>(point);
    int z = get<2>(point);

    bool xNotLow = (x != 0);
    bool yNotLow = (y != 0);
    bool zNotLow = (z != 0);

    bool xNotHigh = (x != sf->size_x-1);
    bool yNotHigh = (y != sf->size_y-1);
    bool zNotHigh = (z != sf->size_z-1);

    if( xNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x-1,y,z), component, components )){
        heap->push({x-1,y,z});
        visitedTracker.setMark( sf->coordsToIndex(x-1,y,z), component, components );
    }

    if( xNotLow && yNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x-1,y+1,z), component, components ) ){
        heap->push({x-1,y+1,z});
        visitedTracker.setMark( sf->coordsToIndex(x-1,y+1,z), component, components );            
    }

    if( yNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x,y+1,z), component, components )){
        heap->push({x,y+1,z});
        visitedTracker.setMark( sf->coordsToIndex(x,y+1,z), component, components );
    }

    if( xNotLow && zNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x-1,y,z+1), component, components )){
        heap->push({x-1,y,z+1});
        visitedTracker.setMark( sf->coordsToIndex(x-1,y,z+1), component, components );
    }

    if( zNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x,y,z+1), component, components )){
        heap->push({x,y,z+1});
        visitedTracker.setMark( sf->coordsToIndex(x,y,z+1), component, components );
    }

    if( yNotHigh && zNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x,y+1,z+1), component, components )){
        heap->push({x,y+1,z+1});
        visitedTracker.setMark( sf->coordsToIndex(x,y+1,z+1), component, components );
    }

    if( xNotLow && yNotHigh && zNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x-1,y+1,z+1), component, components )){
        heap->push({x-1,y+1,z+1});
        visitedTracker.setMark( sf->coordsToIndex(x-1,y+1,z+1), component, components );
    }

    if( xNotHigh && !visitedTracker.checkMark( sf->coordsToIndex(x+1,y,z), component, components )){
        heap->push({x+1,y,z});
        visitedTracker.setMark( sf->coordsToIndex(x+1,y,z), component, components );
    }

    if( xNotHigh && yNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x+1,y-1,z), component, components )){
        heap->push({x+1,y-1,z});
        visitedTracker.setMark( sf->coordsToIndex(x+1,y-1,z), component, components );
    }

    if( yNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x,y-1,z), component, components )){
        heap->push({x,y-1,z});
        visitedTracker.setMark( sf->coordsToIndex(x,y-1,z), component, components );
    }

    if( xNotHigh && zNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x+1,y,z-1), component, components )){
        heap->push({x+1,y,z-1});
        visitedTracker.setMark( sf->coordsToIndex(x+1,y,z-1), component, components );
    }

    if( zNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x,y,z-1), component, components )){
        heap->push({x,y,z-1});
        visitedTracker.setMark( sf->coordsToIndex(x,y,z-1), component, components );
    }

    if( yNotLow && zNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x,y-1,z-1), component, components )){
        heap->push({x,y-1,z-1});
        visitedTracker.setMark( sf->coordsToIndex(x,y-1,z-1), component, components );
    }

    if( xNotHigh && yNotLow && zNotLow && !visitedTracker.checkMark( sf->coordsToIndex(x+1,y-1,z-1), component, components )){
        heap->push({x+1,y-1,z-1});
        visitedTracker.setMark( sf->coordsToIndex(x+1,y-1,z-1), component, components );
    }
}

// the same as before, but does not use an array. This is for when we are just passing over subsets
// and cannot guarantee a traversal order of the subsets.
template<typename T>
inline tuple<bool, bool> isMinOrMax(ScalarField<T>* sf, int i, int j, int k){

    bool isMin = true;
    bool isMax = true;

    if( isMin || isMax ){

        bool iNotLow = (i != 0);
        
        if(iNotLow){
            if(sf->less(i-1,j,k,i,j,k)){
                isMin = false;
            } else {
                isMax = false;
            }
        }

        if( isMin || isMax ){

            bool jNotHigh = (j != sf->size_y-1);

            if(iNotLow && jNotHigh){
                if(sf->less(i-1,j+1,k,i,j,k)){
                    isMin = false;
                } else {
                    isMax = false;
                }
            }

            if( isMin || isMax ){

                if(jNotHigh){
                    if(sf->less(i,j+1,k,i,j,k)){
                        isMin = false;
                    } else {
                        isMax = false;
                    }
                }

                if( isMin || isMax ){

                    bool kNotHigh = (k != sf->size_z-1);

                    if(iNotLow && kNotHigh){
                        if(sf->less(i-1,j,k+1,i,j,k)){
                            isMin = false;
                        } else {
                            isMax = false;
                        }
                    }

                    if( isMin || isMax ){

                        if(kNotHigh){
                            if(sf->less(i,j,k+1,i,j,k)){
                                isMin = false;
                            } else {
                                isMax = false;
                            }
                        }

                        if( isMin || isMax ){

                            if(jNotHigh && kNotHigh){
                                if(sf->less(i,j+1,k+1,i,j,k)){
                                    isMin = false;
                                } else {
                                    isMax = false;
                                }
                            }

                            if( isMin || isMax ){

                                if(iNotLow && jNotHigh && kNotHigh){
                                    if(sf->less(i-1,j+1,k+1,i,j,k)){
                                        isMin = false;
                                    } else {
                                        isMax = false;
                                    }
                                }

                                if( isMin || isMax ){

                                    bool iNotHigh = (i != sf->size_x - 1);

                                    if(iNotHigh){
                                        if(sf->less(i+1,j,k,i,j,k)){
                                            isMin = false;
                                        } else {
                                            isMax = false;
                                        }
                                    }

                                    if( isMin || isMax ){

                                        bool jNotLow = (j != 0);

                                        if(iNotHigh && jNotLow){
                                            if(sf->less(i+1,j-1,k,i,j,k)){
                                                isMin = false;
                                            } else {
                                                isMax = false;
                                            }
                                        }

                                        if( isMin || isMax ){

                                            if(jNotLow){
                                                if(sf->less(i,j-1,k,i,j,k)){
                                                    isMin = false;
                                                } else {
                                                    isMax = false;
                                                }
                                            }

                                            if( isMin || isMax ){

                                                bool kNotLow = (k != 0);

                                                if(iNotHigh && kNotLow){
                                                    if(sf->less(i+1,j,k-1,i,j,k)){
                                                        isMin = false;
                                                    } else {
                                                        isMax = false;
                                                    }
                                                }

                                                if( isMin || isMax ){

                                                    if(kNotLow){
                                                        if(sf->less(i,j,k-1,i,j,k)){
                                                            isMin = false;
                                                        } else {
                                                            isMax = false;
                                                        }
                                                    }

                                                    if( isMin || isMax ){

                                                        if(jNotLow && kNotLow){
                                                            if(sf->less(i,j-1,k-1,i,j,k)){
                                                                isMin = false;
                                                            } else {
                                                                isMax = false;
                                                            }
                                                        }

                                                        if( isMin || isMax ){

                                                            if( iNotHigh && jNotLow && kNotLow ){
                                                                if(sf->less(i+1,j-1,k-1,i,j,k)){
                                                                    isMin = false;
                                                                } else {
                                                                    isMax = false;
                                                                }
                                                            }

                                                        }

                                                    }

                                                }

                                            }

                                        }

                                    }

                                }

                            }

                        }

                    }

                }

            }

        }

    }

    return tuple<bool, bool>({isMin, isMax});
}

// this is really for the region class, but it involves an unrolled loop.
// returns a bit code for the given point, which is a member of the provided set of poitns,
// where each bit corresponds to a neighboring point. That bit is set to 1 if that neighbor is not
// in the given set, and set to 0 if the neighbor is in the set (and thus, the code will be 0)
// if the point is in the interior, and nonzero if it is on the boundary.

// The order is just the order defined here. Not much to add really.
template <typename T>
inline int getBoundaryCode( int point, unordered_set<int>& points, ScalarField<T>* sf ){

    int boundaryCode = 0;

    tuple<int, int, int> coords = sf->indexToCoords(point);

    int yShift = sf->size_x;
    int zShift = sf->size_x*sf->size_y;

    bool xNotLow = ( get<0>(coords) != 0 );
    bool xNotHigh = ( get<0>(coords) != sf->size_x-1 );
    bool yNotLow = ( get<1>(coords) != 0 );
    bool yNotHigh = ( get<1>(coords) != sf->size_y-1 );
    bool zNotLow = ( get<2>(coords) != 0 );
    bool zNotHigh = ( get<2>(coords) != sf->size_z-1 );

    if( xNotLow && !points.count( point - 1 ) ){
        boundaryCode |= (1 << 0);
    }

    if( xNotLow && yNotHigh && !points.count( point + yShift - 1 ) ){
        boundaryCode |= (1 << 1);
    }

    if( yNotHigh && !points.count( point + yShift ) ){
        boundaryCode |= (1 << 2);
    }

    if( xNotLow && zNotHigh && !points.count( point + zShift - 1 ) ){
        boundaryCode |= (1 << 3);
    }

    if( zNotHigh && !points.count( point + zShift ) ){
        boundaryCode |= (1 << 4);
    }

    if( yNotHigh && zNotHigh && !points.count( point + yShift + zShift ) ){
        boundaryCode |= (1 << 5);
    }

    if( xNotLow && yNotHigh && zNotHigh && !points.count( point + yShift + zShift - 1) ){
        boundaryCode |= (1 << 6);
    }

    if( xNotHigh && !points.count( point + 1 ) ){
        boundaryCode |= (1 << 7);
    }

    if( xNotHigh && yNotLow && !points.count( point - yShift + 1 ) ){
        boundaryCode |= (1 << 8);
    }

    if( yNotLow && !points.count( point - yShift ) ){
        boundaryCode |= (1 << 9);
    }

    if( xNotHigh && zNotLow && !points.count( point - zShift + 1 ) ){
        boundaryCode |= (1 << 10);
    }

    if( zNotLow && !points.count( point - zShift ) ){
        boundaryCode |= (1 << 11);
    }

    if( yNotLow && zNotLow && !points.count( point - yShift - zShift ) ){
        boundaryCode |= (1 << 12);
    }

    if( xNotHigh && yNotLow && zNotLow && !points.count( point - yShift - zShift + 1 ) ){
        boundaryCode |= (1 << 13);
    }

    return boundaryCode;

}

template <typename T>
inline void addNeighborsToRegion( int point, int code, unordered_set<int>& pointSet, vector<int>& pointList, ScalarField<T>* sf ){
    
    tuple<int, int, int> coords = sf->indexToCoords(point);

    int yShift = sf->size_x;
    int zShift = sf->size_x*sf->size_y;

    if( (code & (1 << 0)) && !pointSet.count(point - 1) ){
        pointSet.insert(point - 1);
        pointList.push_back(point - 1);
    }

    if( (code & (1 << 1)) && !pointSet.count(point + yShift - 1) ){
        pointSet.insert(point + yShift - 1);
        pointList.push_back(point + yShift - 1);
    }

    if( (code & (1 << 2)) && !pointSet.count(point + yShift) ){
        pointSet.insert(point + yShift);
        pointList.push_back(point + yShift);
    }

    if( (code & (1 << 3)) && !pointSet.count(point + zShift - 1) ){
        pointSet.insert(point + zShift - 1);
        pointList.push_back(point + zShift - 1);
    }

    if( (code & (1 << 4)) && !pointSet.count(point + zShift) ){
        pointSet.insert(point + zShift);
        pointList.push_back(point + zShift);
    }

    if( (code & (1 << 5)) && !pointSet.count(point + yShift + zShift) ){
        pointSet.insert(point + yShift + zShift);
        pointList.push_back(point + yShift + zShift);
    }

    if( (code & (1 << 6)) && !pointSet.count(point + yShift + zShift - 1) ){
        pointSet.insert(point + yShift + zShift - 1);
        pointList.push_back(point + yShift + zShift - 1);
    }

    if( (code & (1 << 7)) && !pointSet.count(point + 1) ){
        pointSet.insert(point + 1);
        pointList.push_back(point + 1);
    }

    if( (code & (1 << 8)) && !pointSet.count(point - yShift + 1) ){
        pointSet.insert(point - yShift + 1);
        pointList.push_back(point - yShift + 1);
    }

    if( (code & (1 << 9)) && !pointSet.count(point - yShift) ){
        pointSet.insert(point - yShift);
        pointList.push_back(point - yShift);
    }

    if( (code & (1 << 10)) && !pointSet.count(point - zShift + 1) ){
        pointSet.insert(point - zShift + 1);
        pointList.push_back(point - zShift + 1);
    }

    if( (code & (1 << 11)) && !pointSet.count(point - zShift) ){
        pointSet.insert(point - zShift);
        pointList.push_back(point - zShift);
    }

    if( (code & (1 << 12)) && !pointSet.count(point - yShift - zShift) ){
        pointSet.insert(point - yShift - zShift);
        pointList.push_back(point - yShift - zShift);
    }

    if( (code & (1 << 13)) && !pointSet.count(point - yShift - zShift + 1) ){
        pointSet.insert(point - yShift - zShift + 1);
        pointList.push_back(point - yShift - zShift + 1);
    }

}