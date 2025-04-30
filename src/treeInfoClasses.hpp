#pragma once

#include <unordered_map>
#include <string>

#include "criticalpointstorage.hpp"
#include "criticalpointtree.hpp"
#include "multitracker.hpp"

template <typename T> 
class MergeTreeInfo{
    public:
        CriticalPointTree<T>* tree;
        vector<shared_ptr<CriticalPoint<T>>>* vertices;
        MultiTracker* visitedTracker;
        unordered_map<int, int>* visitedTrackerIdToSegmentation;        
        vector<int>* USSC; // unvisited saddle segmentation cutoffs
        vector<int>* USCI; // unvisited saddle segmentation indices

        MergeTreeInfo( CriticalPointTree<T>* tree_, vector<shared_ptr<CriticalPoint<T>>>* vertices_, MultiTracker* visitedTracker_, unordered_map<int, int>* visitedTrackerIdToSegmentation_, vector<int>* USSC_, vector<int>* USCI_ ){
            tree = tree_;
            vertices = vertices_;
            visitedTracker = visitedTracker_;
            visitedTrackerIdToSegmentation = visitedTrackerIdToSegmentation_;
            USSC = USSC_;
            USCI = USCI_;
        }

        ~MergeTreeInfo(){
            delete tree;
            delete vertices;
            delete visitedTracker;
            delete visitedTrackerIdToSegmentation;
            delete USSC;
            delete USCI;
        }
};

template <typename T>
class ContourTreeInfo{

    public:
        CriticalPointTree<T>* contourTree = NULL;
        MergeTreeInfo<T>* joinTreeInfo = NULL;
        MergeTreeInfo<T>* splitTreeInfo = NULL;
        vector<int>* contourTreeNodes;
        unordered_map<int, unordered_map<int, pair<T,T>>>* segmentBounds;
        int unsimplifiedSize;

        ContourTreeInfo( CriticalPointTree<T>* contourTree_, vector<int>* contourTreeNodes_, MergeTreeInfo<T>* joinTreeInfo_, MergeTreeInfo<T>* splitTreeInfo_, unordered_map<int, unordered_map<int, pair<T,T>>>* segmentBounds_, int unsimplifiedSize_ ){
            contourTree = contourTree_;
            contourTreeNodes = contourTreeNodes_;
            joinTreeInfo = joinTreeInfo_;
            splitTreeInfo = splitTreeInfo_;
            segmentBounds = segmentBounds_;
            unsimplifiedSize = unsimplifiedSize_;
        }

        void saveToVtu(string filename){
            int numCP = contourTreeNodes->size();

            unordered_map<int, int> globalIndexToLocalIndex;
            int localIndexToGlobalIndex[numCP];

            int i = 0;
            for(auto id : *contourTreeNodes){
                localIndexToGlobalIndex[i] = id;
                globalIndexToLocalIndex[id] = i;
                ++i;
            }

            int x[numCP];
            int y[numCP];
            int z[numCP];

            for(int i = 0; i < numCP; ++i){
                shared_ptr<CriticalPoint<T>> cp = contourTree->getNode( localIndexToGlobalIndex[i] );
                x[i] = cp->x;
                y[i] = cp->y;
                z[i] = cp->z;
            }

            vector<pair<int,int>> edges;
            edges.reserve(numCP - 1);

            for( int i = 0; i < numCP; i++ ){
                int parentGlobalIndex = contourTree->getParentGlobalIndex( localIndexToGlobalIndex[i] );
                if( parentGlobalIndex != -1 ){
                    edges.push_back( {i, globalIndexToLocalIndex[ parentGlobalIndex ]} );
                }
            }

            saveGraph(filename, x, y, z, numCP, edges);            
        }

        ~ContourTreeInfo(){
            delete contourTree;
            delete joinTreeInfo;
            delete splitTreeInfo;
            delete contourTreeNodes;
            delete segmentBounds;
        }

};

template <typename T>
class ContourTreeInfoWithSegmentation{
    public:
        ContourTreeInfo<T>* contourTreeInfo;
        vector<vector<int>>* CTSegmentationIdToPointList;
        unordered_map<int,int>* leafIdToCTSegmentationId;
        vector<pair<T,T>>* bounds;

        ContourTreeInfoWithSegmentation(ContourTreeInfo<T>* contourTreeInfo_, vector<vector<int>>* CTSegmentationIdToPointList_, unordered_map<int,int>* leafIdToCTSegmentationId_, vector<pair<T,T>>* bounds_){
            contourTreeInfo = contourTreeInfo_;
            CTSegmentationIdToPointList = CTSegmentationIdToPointList_;
            leafIdToCTSegmentationId = leafIdToCTSegmentationId_;
            bounds = bounds_;
        }

        ~ContourTreeInfoWithSegmentation(){
            delete contourTreeInfo;
            delete CTSegmentationIdToPointList;
            delete leafIdToCTSegmentationId;
            delete bounds;
        }
};