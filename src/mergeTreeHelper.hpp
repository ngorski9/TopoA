#pragma once

#include <utility>
#include <vector>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>

#include "evtk.hpp"
#include "scalarfield.hpp"
#include "scalarfieldheap.hpp"
#include "unionfind.hpp"
#include "multitracker.hpp"
#include "unrolledloops.hpp"
#include "criticalpointstorage.hpp"
#include "criticalpointtree.hpp"
#include "treeInfoClasses.hpp"

using namespace std;

template <typename T>
pair<MergeTreeInfo<T>*, MergeTreeInfo<T>*> getJoinAndSplitTrees(ScalarField<T>* sf, T epsilon, bool verbose = false){

    int totalPoints = sf->size();

    vector<shared_ptr<CriticalPoint<T>>> mins;
    deque<int> activeMinCriticalPoints;
    int minIndex = 0;

    vector<shared_ptr<CriticalPoint<T>>> maxs;
    deque<int> activeMaxCriticalPoints;
    int maxIndex = 0;

    for( int k = 0; k < sf->size_z; ++k ){
        for( int j = 0; j < sf->size_y; ++j ){
            for( int i = 0; i < sf->size_x; ++i ){
                tuple<bool, bool> minOrMax = isMinOrMax(sf, i, j, k);

                if(get<0>(minOrMax)){
                    mins.push_back( shared_ptr<CriticalPoint<T>>(new CriticalPoint<T>(true, minIndex, i, j, k, sf)) );
                    activeMinCriticalPoints.push_back(minIndex);
                    minIndex += 1;
                } else if(get<1>(minOrMax)){
                    maxs.push_back( shared_ptr<CriticalPoint<T>>(new CriticalPoint<T>(false, maxIndex, i, j, k, sf)) );
                    activeMaxCriticalPoints.push_back(maxIndex);
                    maxIndex += 1;
                }

            } // end iter x
        } // end iter y
    } // end iter zsh

    // needed later so that we can distinguish between leaf nodes and interior nodes
    int numMinLeafs = mins.size();
    int numMaxLeafs = maxs.size();

    if( verbose ){
        cout << "found extrema" << endl;
    }

    // Join tree & Split tree

    MultiTracker* minSegmentation = new MultiTracker(sf->size());
    MultiTracker* maxSegmentation = new MultiTracker(sf->size());

    vector<shared_ptr<CriticalPoint<T>>> unvisitedMinSaddles;
    vector<shared_ptr<CriticalPoint<T>>> unvisitedMaxSaddles;

    CriticalPointLess<T> minComparison;
    CriticalPointMore<T> maxComparison;

    unordered_map<int, int> globalIndexToMinIndex;
    unordered_map<int, int> globalIndexToMaxIndex;

    minComparison.setSF(sf);
    maxComparison.setSF(sf);

    // Critical type 0 = join tree. Critical type 1 = split tree.
    for( int criticalType = 0; criticalType < 2; criticalType++ ){
        // variables whose value depends on if we are doing a min or max
        bool min = (criticalType == 0);
        vector<shared_ptr<CriticalPoint<T>>>& cps = min ? mins : maxs;
        MultiTracker* visitedTracker = min ? minSegmentation : maxSegmentation;
        deque<int>& activeCriticalPoints = min ? activeMinCriticalPoints : activeMaxCriticalPoints;
        vector<shared_ptr<CriticalPoint<T>>>& unvisitedSaddles = min ? unvisitedMinSaddles : unvisitedMaxSaddles;
        unordered_map<int, int>& globalIndexToCriticalIndex = min ? globalIndexToMinIndex : globalIndexToMaxIndex;

        int numCritical = cps.size();
        int numComponents = numCritical;
        int maxNumCritical = 2*numCritical-1;

        // add extrema to the globalIndexToCriticalIndex
        for( int i = 0; i < numCritical; i++ ){
            globalIndexToCriticalIndex[ cps[i]->globalIndex ] = i;
        }

        cps.reserve(maxNumCritical);

        MultiTracker heapTracker = MultiTracker(sf->size());
        UnionFind components = UnionFind(maxNumCritical);

        while( numComponents > 1 && activeCriticalPoints.size() > 1){
            int cpIndex = activeCriticalPoints.front();

            activeCriticalPoints.pop_front();
            shared_ptr<CriticalPoint<T>> cp = cps[cpIndex];
            if( cp->criticalType == 0 || cp->criticalType == 3 ){
                heapTracker.setMark( sf->coordsToIndex( cp->x, cp->y, cp->z ), cpIndex, components );
            }

            tuple<int,int,int> lastTuple = cp->heap->front();

            bool keepGoing = true;
            while( keepGoing ){
                if(cp->heap->size == 0){
                    keepGoing = false;
                    cout << "empty heap" << endl;
                    break;
                }
                tuple<int,int,int> nextTuple = cp->heap->pop();

                // sadly wee need to double up the check here because equality breaks it.
                if( (min && sf->less( get<0>(nextTuple), get<1>(nextTuple), get<2>(nextTuple), get<0>(lastTuple), get<1>(lastTuple), get<2>(lastTuple) ))
                || (!min && sf->more( get<0>(nextTuple), get<1>(nextTuple), get<2>(nextTuple), get<0>(lastTuple), get<1>(lastTuple), get<2>(lastTuple) ))){
                    // last tuple is a saddle point
                    keepGoing = false;

                    int newSaddleGlobalIndex = sf->coordsToIndex( get<0>(lastTuple), get<1>(lastTuple), get<2>(lastTuple) );

                    if( globalIndexToCriticalIndex.find( newSaddleGlobalIndex ) != globalIndexToCriticalIndex.end() ){
                        
                        numComponents -= 1;
                        // this is not the first time the saddle has been discovered
                        int newCriticalIndex = globalIndexToCriticalIndex[ newSaddleGlobalIndex ];
                        components.union_(cpIndex, newCriticalIndex);
                        shared_ptr<CriticalPoint<T>> newCP = cps[newCriticalIndex];
                        newCP->children.push_back(cp);

                        cp->parent = newCP;
                        newCP->mergeInHeap(cp);
                                                                    // dereference here to not break the old  code.
                        if( !newCP->active && newCP->linkIsVisited(*visitedTracker, components) ){
                            newCP->active = true;
                            activeCriticalPoints.push_back(newCriticalIndex);
                        }

                    } else {
                        // this is the first time hte saddle has been discovered
                        globalIndexToCriticalIndex[ sf->coordsToIndex( get<0>(lastTuple), get<1>(lastTuple), get<2>(lastTuple) ) ] = numCritical;
                        shared_ptr<CriticalPoint<T>> newCP = shared_ptr<CriticalPoint<T>>(new CriticalPoint<T>(min, numCritical, get<0>(lastTuple), get<1>(lastTuple), get<2>(lastTuple), cp->heap, sf));
                        newCP->children.push_back(cp);
                        cps.push_back(newCP);
                        components.union_(cpIndex, numCritical);
                        cp->heap = NULL;
                        cp->parent = newCP;
                        numCritical += 1;
                    }

                } else {
                    // otherwise, add all neighboring points to the queue
                    // we leave  this for the loop unrolling file.
                    visitedTracker->setMark( sf->coordsToIndex( get<0>(nextTuple), get<1>(nextTuple), get<2>(nextTuple) ), cpIndex, components );

                    lastTuple = nextTuple;
                    addNeighborsToHeap(nextTuple, cp->heap, sf, heapTracker, components, cpIndex);
                }

            }

        }

        for(int i = 0; i < numCritical; ++i){
            if( !cps[i]->active ){
                unvisitedSaddles.push_back( cps[i] );
            }
        }

        for( const int& idx : activeCriticalPoints ){
            unvisitedSaddles.push_back( cps[idx] );
        }

        if( criticalType == 0 ){
            sort( unvisitedSaddles.begin(), unvisitedSaddles.end(), minComparison );
        } else {
            sort( unvisitedSaddles.begin(), unvisitedSaddles.end(), maxComparison );
        }

        int numUnvisited = unvisitedSaddles.size();
        for(int i = 1; i < numUnvisited; ++i){
            unvisitedSaddles[i-1]->parent = unvisitedSaddles[i];
            unvisitedSaddles[i]->children.push_back(unvisitedSaddles[i-1]);
        }

        if( verbose ){
            if( criticalType == 0 ){
                cout << "computed join tree" << endl;
            } else {
                cout << "computed split tree" << endl;
            }
        }

    }

    // define and simplify the join and split tree

    CriticalPointTree<T>* joinTree = new CriticalPointTree<T>(sf);
    CriticalPointTree<T>* splitTree = new CriticalPointTree<T>(sf);

    vector<shared_ptr<CriticalPoint<T>>> processStack;
    // critical type 0 = join tree, 1 = split tree
    for( int criticalType = 0; criticalType < 2; criticalType++ ){

        bool isJoinTree = (criticalType == 0);
        vector<shared_ptr<CriticalPoint<T>>>& cpList = isJoinTree ? mins : maxs;
        CriticalPointTree<T>& tree = isJoinTree ? (*joinTree) : (*splitTree);

        for( auto m : cpList ){

            while( m != NULL && !tree.hasNode( m->globalIndex ) ){
                processStack.push_back(m);
                m = m->parent;
            }

            while( processStack.size() > 0 ){
                shared_ptr<CriticalPoint<T>>& nextNode = processStack.back();

                tree.newNode( nextNode->criticalType, nextNode->criticalIndex, nextNode->x, nextNode->y, nextNode->z );
                if( nextNode->parent != NULL ){
                    tree.connectNodes( nextNode->globalIndex, nextNode->parent->globalIndex );
                }

                processStack.pop_back();
            }
        }
    }

    // simplify these trees

    int numMins = mins.size();
    int numMaxs = maxs.size();

    vector<shared_ptr<CriticalPoint<T>>>* minsSorted = new vector<shared_ptr<CriticalPoint<T>>>();
    vector<shared_ptr<CriticalPoint<T>>>* maxsSorted = new vector<shared_ptr<CriticalPoint<T>>>();

    minsSorted->reserve(numMins);
    maxsSorted->reserve(numMaxs);

    for( auto m : mins ){
        minsSorted->push_back( joinTree->getNode( m->globalIndex ) );
    }

    for( auto m : maxs ){
        maxsSorted->push_back( splitTree->getNode( m->globalIndex ) );
    }

    // sort mins from greatest to smallest, maxes from smallest to greatest.
    sort( minsSorted->begin(), minsSorted->begin() + numMinLeafs, maxComparison );
    sort( maxsSorted->begin(), maxsSorted->begin() + numMaxLeafs, minComparison );

    // simplify the tree and track how the segmentation merges.
    UnionFind jtMergedSegmentation(numMins);
    UnionFind stMergedSegmentation(numMaxs);

    for( int criticalType = 0; criticalType < 2; criticalType += 1 ){

        bool isJT = (criticalType == 0);

        int& numLeafs = isJT ? numMinLeafs : numMaxLeafs;
        vector<shared_ptr<CriticalPoint<T>>>& leafs = isJT ? (*minsSorted) : (*maxsSorted);
        UnionFind& mergedSegmentation = isJT ? jtMergedSegmentation : stMergedSegmentation;
        CriticalPointTree<T>& tree = isJT ? (*joinTree) : (*splitTree);

        int i = 0;

        while(i < numLeafs - 1){

            int prunedLeafGlobalIndex = leafs[i]->globalIndex;
            if(tree.getEdgeLength( prunedLeafGlobalIndex ) < epsilon){

                int prunedLeafSegmentationIndex  = leafs[i]->criticalIndex;

                shared_ptr<CriticalPoint<T>> parent = leafs[i]->parent;
                tree.deleteNode( prunedLeafGlobalIndex );
                leafs.erase(leafs.begin() + i);

                mergedSegmentation.union_( prunedLeafSegmentationIndex, parent->children[0]->criticalIndex );

                if( parent->children.size() == 1 ){
                    mergedSegmentation.union_( parent->criticalIndex, parent->children[0]->criticalIndex );

                    leafs.erase( find(leafs.begin(), leafs.end(), parent) );
                    tree.deleteNode( parent->globalIndex );
                }

                --numLeafs;
            } else {
                ++i;
            }

        }
    }

    // order the unvisited saddles such that indexing into the list shows the segmentation in the unvisited tree.
    // this prevents us from needing to compute the segmentation; we can just compute it when the error bounds are computed later.

    vector<int>* jtUSSC = new vector<int>(); // join tree unvisited saddle segmentation cutoffs
    vector<int>* jtUSCI = new vector<int>(); // join tree unvisited saddle critical indices
    vector<int>* stUSSC = new vector<int>(); // split tree unvisited saddle segmentation cutoffs
    vector<int>* stUSCI = new vector<int>(); // split tree unvisited saddle critical indices

    for( int criticalType = 0; criticalType < 2; criticalType++ ){

        bool isJT = (criticalType == 0);

        vector<int>& USSC = (isJT) ? (*jtUSSC) : (*stUSSC);
        vector<int>& USCI = (isJT) ? (*jtUSCI) : (*stUSCI);
        vector<shared_ptr<CriticalPoint<T>>>& unmergedSaddles = (isJT) ? unvisitedMinSaddles : unvisitedMaxSaddles;
        UnionFind& mergedSegmentation = (isJT) ? jtMergedSegmentation : stMergedSegmentation;

        int lastSegmentation = -1;
        int numUnmergedSaddles = unmergedSaddles.size();

        USSC.reserve(numUnmergedSaddles);
        USCI.reserve(numUnmergedSaddles);

        for(int i = 0; i < numUnmergedSaddles; ++i){

            int segmentation = mergedSegmentation.find_( unmergedSaddles[i]->criticalIndex );

            if( segmentation != lastSegmentation ){
                USSC.push_back( unmergedSaddles[i]->globalIndex );
                USCI.push_back( unmergedSaddles[i]->criticalIndex );
                lastSegmentation = segmentation;
            }

        }

    }

    unordered_map<int, int>* jtVTIS = new unordered_map<int,int>(); // join tree visited tracker id to segmentation
    unordered_map<int, int>* stVTIS = new unordered_map<int,int>(); // split tree visited tracker id to segmentation

    for( int criticalType = 0; criticalType < 2; ++criticalType ){

        bool isJT = (criticalType == 0);

        unordered_map<int, int>& VTIS = isJT ? (*jtVTIS) : (*stVTIS);
        vector<shared_ptr<CriticalPoint<T>>>& treePoints = isJT ? (*minsSorted) : (*maxsSorted);
        vector<shared_ptr<CriticalPoint<T>>>& allMTPoints = isJT ? mins : maxs;
        UnionFind& mergedSegmentation = isJT ? jtMergedSegmentation : stMergedSegmentation;

        int numMTPoints = allMTPoints.size();
        for( int i = 0; i < numMTPoints; ++i ){
            VTIS[ i ] = allMTPoints[ mergedSegmentation.find_( i ) ]->globalIndex;
        }

    }

    // delete heaps as they will not be needed any more.
    for( auto m : (*minsSorted) ){
        delete m->heap;
    }

    for( auto m : (*maxsSorted) ){
        delete m->heap;
    }

    MergeTreeInfo<T>* joinTreeInfo = new MergeTreeInfo<T>( joinTree, minsSorted, minSegmentation, jtVTIS, jtUSSC, jtUSCI );
    MergeTreeInfo<T>* splitTreeInfo = new MergeTreeInfo<T>( splitTree, maxsSorted, maxSegmentation, stVTIS, stUSSC, stUSCI );

    return { joinTreeInfo, splitTreeInfo };

} // end function