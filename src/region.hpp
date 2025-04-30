#pragma once

#include <unordered_set>
#include <vector>
#include <memory>

#include "unrolledloops.hpp"
#include "scalarfield.hpp"
#include "compressorData.hpp"

using namespace std; // :skull:

template<typename T>
class Region{

    private:
    public:
        vector<int> boundaryPoints;
        vector<int> boundaryCodes;
        
        unordered_set<int> pointSet;

        ScalarField<T>* sf;

        void setupBoundary(){

            boundaryPoints.reserve( pointList.size() );
            boundaryCodes.reserve( pointList.size() );

            // calculate boundary
            for( auto i : pointList ){
                
                int boundaryCode = getBoundaryCode<T>( i, pointSet, sf );
                if( boundaryCode != 0 ){
                    boundaryPoints.push_back(i);
                    boundaryCodes.push_back(boundaryCode);
                }

            }
        }

    public:
        vector<int> pointList;

        Region(int point, ScalarField<T>* sf_ ){
            sf = sf_;
            pointList = {point};
            pointSet.insert(point);
            boundaryPoints = {point};
            boundaryCodes = {getBoundaryCode<T>(point, pointSet, sf)};
        }

        Region( vector<int> pointList_, ScalarField<T>* sf_ ){
            pointList.reserve(pointList_.size());

            for( auto i : pointList_ ){
                if( pointSet.count(i) == 0 ){
                    pointSet.insert(i);
                    pointList.push_back(i);
                }
            }

            sf = sf_;

            setupBoundary();
        }

        Region( vector<shared_ptr<vector<int>>> pointListList, ScalarField<T>* sf_ ){
            sf = sf_;
            int total = 0;

            for( auto l : pointListList ){
                total += l->size();
            }
            pointList.reserve(total);

            for( auto l : pointListList ){
                for( auto i : (*l) ){
                    if(pointSet.count(i) == 0){
                        pointSet.insert(i);
                        pointList.push_back(i);
                    }
                }
            }

            setupBoundary();
        }

        bool contains(int idx){
            return pointSet.count(idx) == 1;
        }

        // grows the region by expanding around a single layer
        void grow(){

            int originalSize = pointList.size();
            int boundarySize = boundaryPoints.size();
            for( int i = 0; i < boundarySize; ++i){
                addNeighborsToRegion<T>( boundaryPoints[i], boundaryCodes[i], pointSet, pointList, sf );
            }

            boundaryPoints.clear();
            boundaryCodes.clear();

            // the boundary points of the new region will simply be the most recently added points.
            // also use the newly added points to update the min and max of the region.
            int newSize = pointList.size();
            for(int i = originalSize; i < newSize; ++i){
                int boundaryCode = getBoundaryCode<T>( pointList[i], pointSet, sf );
                if( boundaryCode != 0 ){
                    boundaryPoints.push_back(pointList[i]);
                    boundaryCodes.push_back(boundaryCode);
                }

            }

        }

        // just repeats the grow function n times.
        void grow(int n){

            for(int i = 0; i < n; ++i){
                grow();
            }

        }

        void mergeInOtherRegion(Region<T>& other){

            // combine the other region into the current one.
            // do not do anything fancy with the boundary here. Just combine them.
            // This isn't the exact boundary, but it will correct itself after one layer of growth.
            // and the intended use case is to merge in a small ball around saddle points.

            // also we do not reserve here because we intend for "other" to be smaller, in which case
            // the size of both point vectors combined is less than double the size of the current memory vector.
            // we would like for the memory of the current one to double, if needed, because we will be expanding the list as we go.

            // because boundary points are added in the same order that points are added,
            // we check if an added point is a boundary point by checking if each individual point is equal
            // to the next boundary point that we have not checked yet.
            int otherNextBoundaryPoint = 0;
            int otherBoundarySize = other.boundaryPoints.size();

            for( int i = 0; i < other.pointList.size(); ++i ){
                if( pointSet.count( other.pointList[i] ) == 0 ){

                    pointSet.insert( other.pointList[i] );
                    pointList.push_back( other.pointList[i] );

                }

                if( otherNextBoundaryPoint < otherBoundarySize && other.boundaryPoints[otherNextBoundaryPoint] == other.pointList[i] ){
                    boundaryPoints.push_back(other.boundaryPoints[otherNextBoundaryPoint]);
                    boundaryCodes.push_back(other.boundaryCodes[otherNextBoundaryPoint]);
                    ++otherNextBoundaryPoint;
                }                
            }
        }

};