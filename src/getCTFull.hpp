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
ContourTreeInfo<T> getContourTree(ScalarField<T>* sf, double epsilon_, bool verbose=false){

    int totalPoints = sf->size();

    vector<shared_ptr<CriticalPoint<T>>> mins;
    deque<int> activeMinCriticalPoints;
    int minIndex = 0;

    vector<shared_ptr<CriticalPoint<T>>> maxs;
    deque<int> activeMaxCriticalPoints;
    int maxIndex = 0;

    T globalMax = -INFINITY;
    T globalMin = INFINITY;
    int globalMaxIndex; // id number of the point where the global max is achieved
    int globalMinIndex; // same but for min

    for( int k = 0; k < sf->size_z; ++k ){
        for( int j = 0; j < sf->size_y; ++j ){
            for( int i = 0; i < sf->size_x; ++i ){
                tuple<bool, bool> minOrMax = isMinOrMax(sf, i, j, k);

                if(get<0>(minOrMax)){
                    mins.push_back( shared_ptr<CriticalPoint<T>>(new CriticalPoint<T>(true, minIndex, i, j, k, sf)) );
                    activeMinCriticalPoints.push_back(minIndex);
                    minIndex += 1;
                    
                    T minValue = sf->getElement(i,j,k);
                    if( minValue < globalMin ){
                        globalMin = minValue;
                        globalMinIndex = sf->coordsToIndex(i,j,k);
                    }
                } else if(get<1>(minOrMax)){
                    maxs.push_back( shared_ptr<CriticalPoint<T>>(new CriticalPoint<T>(false, maxIndex, i, j, k, sf)) );
                    activeMaxCriticalPoints.push_back(maxIndex);
                    maxIndex += 1;

                    T maxValue = sf->getElement(i,j,k);
                    if( maxValue > globalMax ){
                        globalMax = maxValue;
                        globalMaxIndex = sf->coordsToIndex(i,j,k);
                    }
                }

            } // end iter x
        } // end iter y
    } // end iter zsh

    T globalRange = globalMax - globalMin;

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
        MultiTracker& visitedTracker = min ? (*minSegmentation) : (*maxSegmentation);
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

                        if( !newCP->active && newCP->linkIsVisited(visitedTracker, components) ){
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
                    visitedTracker.setMark( sf->coordsToIndex( get<0>(nextTuple), get<1>(nextTuple), get<2>(nextTuple) ), cpIndex, components );

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

    CriticalPointTree<T> joinTree(sf);
    CriticalPointTree<T> splitTree(sf);

    vector<shared_ptr<CriticalPoint<T>>> processStack;
    // critical type 0 = join tree, 1 = split tree
    for( int criticalType = 0; criticalType < 2; criticalType++ ){

        bool isJoinTree = (criticalType == 0);
        vector<shared_ptr<CriticalPoint<T>>>& cpList = isJoinTree ? mins : maxs;
        CriticalPointTree<T>& tree = isJoinTree ? joinTree : splitTree;

        for( auto m : cpList ){

            while( m != NULL && !tree.hasNode( m->globalIndex ) ){
                processStack.push_back(m);
                m = m->parent;
            }

            while( processStack.size() > 0 ){
                shared_ptr<CriticalPoint<T>> nextNode = processStack.back();

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

    int unsimplifiedSize = numMins + numMaxs;

    vector<shared_ptr<CriticalPoint<T>>>* minsSorted = new vector<shared_ptr<CriticalPoint<T>>>();
    vector<shared_ptr<CriticalPoint<T>>>* maxsSorted = new vector<shared_ptr<CriticalPoint<T>>>();

    for( auto m : mins ){
        minsSorted->push_back( joinTree.getNode( m->globalIndex ) );
    }

    for( auto m : maxs ){
        maxsSorted->push_back( splitTree.getNode( m->globalIndex ) );
    }

    // sort mins from greatest to smallest, maxes from smallest to greatest.
    sort( minsSorted->begin(), minsSorted->begin() + numMinLeafs, maxComparison );
    sort( maxsSorted->begin(), maxsSorted->begin() + numMaxLeafs, minComparison );

    bool* keepMins = new bool[mins.size()];
    bool* keepMaxs = new bool[maxs.size()];
    fill_n(keepMins, mins.size(), true);
    fill_n(keepMaxs, maxs.size(), true);
    int numKeptMins = mins.size();
    int numKeptMaxs = maxs.size();

    UnionFind jtMergedSegmentation(numMins);
    UnionFind stMergedSegmentation(numMaxs);

    T epsilon = epsilon_ * globalRange;

    for( int criticalType = 0; criticalType < 2; criticalType += 1 ){

        bool isJT = (criticalType == 0);

        int& numLeafs = isJT ? numMinLeafs : numMaxLeafs;
        vector<shared_ptr<CriticalPoint<T>>>& leafs = isJT ? (*minsSorted) : (*maxsSorted);
        UnionFind& mergedSegmentation = isJT ? jtMergedSegmentation : stMergedSegmentation;
        CriticalPointTree<T>& tree = isJT ? joinTree : splitTree;
        bool* keepLeafs = isJT ? keepMins : keepMaxs;
        int& numLeavesKept = isJT ? numKeptMins : numKeptMaxs;        

        int i = 0;

        while(i < numLeafs){
            int prunedLeafGlobalIndex = leafs[i]->globalIndex;

            if(keepLeafs[leafs[i]->criticalIndex] && tree.getParentGlobalIndex(prunedLeafGlobalIndex) != -1 && tree.getEdgeLength( prunedLeafGlobalIndex ) < epsilon){
                int prunedLeafSegmentationIndex  = leafs[i]->criticalIndex;

                shared_ptr<CriticalPoint<T>> parent = leafs[i]->parent;
                tree.deleteNode( prunedLeafGlobalIndex );
                keepLeafs[leafs[i]->criticalIndex] = false;
                --numLeavesKept;

                mergedSegmentation.union_( prunedLeafSegmentationIndex, parent->children[0]->criticalIndex );

                if( parent->children.size() == 1 ){
                    mergedSegmentation.union_( parent->criticalIndex, parent->children[0]->criticalIndex );

                    keepLeafs[ parent->criticalIndex ] = false;
                    --numLeavesKept;
                    tree.deleteNode( parent->globalIndex );
                }
            }

            ++i;

        }
    }

    if( verbose ){
        cout << "simplified trees" << endl;
    }

    // reset the mins sorted and maxs sorted based on which leaves are kept
    delete minsSorted;
    delete maxsSorted;

    minsSorted = new vector<shared_ptr<CriticalPoint<T>>>();
    maxsSorted = new vector<shared_ptr<CriticalPoint<T>>>();

    minsSorted->reserve(numMinLeafs);
    maxsSorted->reserve(numMaxLeafs);

    for( int i = 0; i < numMins; ++i ){
        if( keepMins[i] ){
            minsSorted->push_back(joinTree.getNode(mins[i]->globalIndex));
        }
    }

    for( int i = 0; i < numMaxs; ++i ){
        if( keepMaxs[i] ){
            maxsSorted->push_back(splitTree.getNode(maxs[i]->globalIndex));
        }
    }

    delete [] keepMins;
    delete [] keepMaxs;

    // clone the join and split trees

    CriticalPointTree<T>* rawJoinTree = new CriticalPointTree<T>(sf);
    CriticalPointTree<T>* rawSplitTree = new CriticalPointTree<T>(sf);

    for( auto node : (*minsSorted) ){
        rawJoinTree->newNode( node->criticalType, node->globalIndex, node->x, node->y, node->z );
    }

    for( auto node : (*minsSorted) ){
        if( node->parent ){
            rawJoinTree->connectNodes( node->globalIndex, node->parent->globalIndex );
        }
    }

    for( auto node : (*maxsSorted) ){
        rawSplitTree->newNode( node->criticalType, node->globalIndex, node->x, node->y, node->z );
    }
    for( auto node : (*maxsSorted) ){
        if( node->parent ){
            rawSplitTree->connectNodes( node->globalIndex, node->parent->globalIndex );
        }
    }

    // create augmented join tree and augmented split tree
    // sort all mins now (not just the root nodes)

    sort( minsSorted->begin(), minsSorted->end(), minComparison );
    sort( maxsSorted->begin(), maxsSorted->end(), maxComparison );

    // exit(1);
    // calculate the segmentation:

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

    if( verbose ){
        cout << "primed segmentation trees" << endl;
    }

    // insert nodes from the join tree into the split tree and split tree into contour tree.
    // the code is very hard to read. I apologize for this.

    // number of critical points that appear in both the jt and the st
    // this can happen in the case of a monkey saddle (otherwise, not possible in 3D)
    int numDoubleCriticalPoints = 0; 

    // keeps track of the segmentation number of each point added to each augmented tree
    unordered_map<int,vector<int>> augmentedJTSegmentation(numMins);
    unordered_map<int,vector<int>> augmentedSTSegmentation(numMaxs);

    ScalarFieldLess<T> scalarFieldLess;
    ScalarFieldMore<T> scalarFieldMore;
    scalarFieldLess.setSF(sf);
    scalarFieldMore.setSF(sf);

    // criticalType = 0: merge join tree into split tree.
    // criticalType = 1: merge split tree into join tree.
    for( int criticalType = 0; criticalType < 2; ++criticalType ){
        bool mergeInJT = (criticalType == 0);

        vector<shared_ptr<CriticalPoint<T>>>& ungrownSaddles = mergeInJT ? unvisitedMaxSaddles : unvisitedMinSaddles;
        vector<shared_ptr<CriticalPoint<T>>>& mergeInPoints = mergeInJT ? (*minsSorted) : (*maxsSorted);
        vector<shared_ptr<CriticalPoint<T>>>& mergeToPoints = mergeInJT ? maxs : mins;
        MultiTracker& segmentationTracker = mergeInJT ? (*maxSegmentation) : (*minSegmentation);
        UnionFind& segmentationMergeTracker = mergeInJT ? stMergedSegmentation : jtMergedSegmentation;
        unordered_map<int,int>& globalIndexToMergingToIndex = mergeInJT ? globalIndexToMaxIndex : globalIndexToMinIndex;
        CriticalPointTree<T>& tree = mergeInJT ? splitTree : joinTree;
        unordered_map<int,vector<int>>& mergedInNewSegmentation = mergeInJT ? augmentedSTSegmentation : augmentedJTSegmentation;
        unordered_map<int,vector<int>>& mergedToNewSegmentation = mergeInJT ? augmentedJTSegmentation : augmentedSTSegmentation;


        vector<int>& USSC = mergeInJT ? (*stUSSC) : (*jtUSSC);
        vector<int>& USCI = mergeInJT ? (*stUSCI) : *(jtUSCI);
        unordered_map<int,int> VTIS = mergeInJT ? (*stVTIS) : (*jtVTIS);

        // iterate through the global indices of points that we are inserting into the tree (to make it augmented)
        for( auto cp : mergeInPoints ){

            // find the node that the inserted node will become its parent by finding the segmentation id of the point we are inserting.
            int segmentation = segmentationTracker.entries[cp->globalIndex];
            int nextChildGlobalIndex; // the global index of the node that we will be inserting above.

            if( segmentation == -1 ){

                if( criticalType == 0 ){
                    nextChildGlobalIndex = VTIS[ USCI[upper_bound( USSC.begin(), USSC.end(), cp->globalIndex, scalarFieldMore ) - USSC.begin() - 1] ];
                } else {
                    nextChildGlobalIndex = VTIS[ USCI[upper_bound( USSC.begin(), USSC.end(), cp->globalIndex, scalarFieldLess ) - USSC.begin() - 1] ];
                }

            } else if( segmentation ==  -2 ) {
                // if the point has multiple segmentation values, break the tie arbitrarily.
                nextChildGlobalIndex = VTIS[segmentationTracker.spillOverMarks[cp->globalIndex][0]];
            } else {
                nextChildGlobalIndex = VTIS[segmentation];
            }

            // if( criticalType == 0 ){
            //     cout << cp->globalIndex << " " << nextChildGlobalIndex << endl;
            // } else {
            //     cout << nextChildGlobalIndex << " " << cp->globalIndex << endl;
            // }

            mergedInNewSegmentation[ cp->globalIndex ] = vector<int>();
            mergedInNewSegmentation[ cp->globalIndex ].push_back(nextChildGlobalIndex);

            mergedToNewSegmentation[ cp->globalIndex ] = vector<int>();
            mergedToNewSegmentation[ cp->globalIndex ].push_back(cp->globalIndex);

            // add the node to the tree above "next child"
            if( !tree.hasNode( cp->globalIndex ) ){
                tree.newNode(cp->criticalType, cp->criticalIndex, cp->x, cp->y, cp->z, nextChildGlobalIndex);
            } else {
                numDoubleCriticalPoints += 1;
            }

        }

    }

    if( verbose ){
        cout << "augmented trees" << endl;
    }

    // divide this by 2 since they are each counted twice.
    numDoubleCriticalPoints /= 2;

    // now that we have the augmented trees, run the algorithm to combine them.
    // the main reason why we do this for the compressor is actually to compute how the segmentation combines for everything.

    unordered_map<int,unordered_map<int,pair<T,T>>>* contourTreeSegmentationBounds = new unordered_map<int,unordered_map<int,pair<T,T>>>();

    vector<int> leafQueue;
    numMins = minsSorted->size();
    numMaxs = maxsSorted->size();
    int numCP = numMins + numMaxs - numDoubleCriticalPoints;
    leafQueue.reserve(numCP);
    int leafQueueHead = 0;
    unordered_set<int> queueMembers;

    for( auto m : (*minsSorted) ){
        if( (joinTree.getNode( m->globalIndex )->children.size() + splitTree.getNode( m->globalIndex )->children.size() == 1) && queueMembers.find(m->globalIndex) == queueMembers.end() ){
            leafQueue.push_back( m->globalIndex );
            queueMembers.insert( m->globalIndex );
        }
    }

    for( auto m : (*maxsSorted) ){
        if( (joinTree.getNode( m->globalIndex )->children.size() + splitTree.getNode( m->globalIndex )->children.size() == 1) && queueMembers.find(m->globalIndex) == queueMembers.end()){
            leafQueue.push_back( m->globalIndex );
            queueMembers.insert( m->globalIndex );
        }
    }

    CriticalPointTree<T>* contourTree = new CriticalPointTree<T>(sf);
    vector<int>* contourTreeNodeIndices = new vector<int>();
    contourTreeNodeIndices->reserve(numCP);

    // used in the loop, but we define these here so that the memory doesn't get repeatedly deleted
    // will store which segments from the join tree and the split tree contribute to the current segment.
    // note that it is ( union of jt segments ) intersect ( union of st segments )
    vector<int> jtContributingSegments;
    vector<int> stContributingSegments;

    int numCombo = 0;

    while( leafQueueHead != numCP-1){
        auto t1 = std::chrono::high_resolution_clock::now();

        int nextLeaf = leafQueue[leafQueueHead];
        ++leafQueueHead;

        // get the branch associated with the leaf that we are adding
        int jtParent = joinTree.getParentGlobalIndex(nextLeaf);
        int stParent = splitTree.getParentGlobalIndex(nextLeaf);
        int nextParent = joinTree.isLeaf(nextLeaf) ? jtParent : stParent;

        if( !contourTree->hasNode( nextLeaf ) ){
            contourTreeNodeIndices->push_back( nextLeaf );
            shared_ptr<CriticalPoint<T>> m = joinTree.getNode( nextLeaf );
            contourTree->newNode( m->criticalType, m->criticalIndex, m->x, m->y, m->z );            
        }
        if( !contourTree->hasNode( nextParent ) ){
            contourTreeNodeIndices->push_back( nextParent );
            shared_ptr<CriticalPoint<T>> m = joinTree.getNode( nextParent );
            contourTree->newNode( m->criticalType, m->criticalIndex, m->x, m->y, m->z );
        }

        // calculate contour tree segmentation bounds
        int smallerNode;
        int largerNode;
        if( scalarFieldLess(nextLeaf, nextParent) ){
            smallerNode = nextLeaf;
            largerNode = nextParent;
        } else {
            smallerNode = nextParent;
            largerNode = nextLeaf;
        }

        jtContributingSegments.clear();
        stContributingSegments.clear();
        unordered_set<int> jtContributingSegmentsSet;
        unordered_set<int> stContributingSegmentsSet;

        shared_ptr<CriticalPoint<T>> smallerNodePtr = joinTree.getNode(smallerNode);
        while( smallerNodePtr != nullptr && smallerNodePtr->globalIndex != largerNode ){
            for( auto segment : augmentedJTSegmentation[smallerNodePtr->globalIndex] ){
                if( !jtContributingSegmentsSet.count(segment) ){
                    jtContributingSegments.push_back(segment);
                    jtContributingSegmentsSet.insert(segment);
                }
            }
            smallerNodePtr = smallerNodePtr->parent;
        }

        shared_ptr<CriticalPoint<T>> largerNodePtr = splitTree.getNode(largerNode);
        while( largerNodePtr != nullptr && largerNodePtr->globalIndex != smallerNode ){
            for( auto segment : augmentedSTSegmentation[largerNodePtr->globalIndex] ){
                if( !stContributingSegmentsSet.count(segment) ){
                    stContributingSegments.push_back(segment);
                    stContributingSegmentsSet.insert(segment);
                }
            }
            largerNodePtr = largerNodePtr->parent;
        }

        T jtUpperBound = -INFINITY;
        T jtLowerBound = INFINITY;
        T stUpperBound = -INFINITY;
        T stLowerBound = INFINITY;

        // based on the contributing segments, calculate the bounds for each combination of join tree or split tree.
        for( auto segment : jtContributingSegments ){
            jtLowerBound = min(sf->getElement(segment), jtLowerBound);

            if( rawJoinTree->getParentGlobalIndex(segment) == -1 ){
                jtUpperBound = INFINITY;
            } else{
                jtUpperBound = max( sf->getElement( rawJoinTree->getParentGlobalIndex(segment) ), jtUpperBound );
            }
        }

        for( auto segment : stContributingSegments ){
            stUpperBound = max( sf->getElement(segment), stUpperBound );
            if( rawSplitTree->getParentGlobalIndex(segment) == -1 ){
                stLowerBound = -INFINITY;
            } else {
                stLowerBound = min( sf->getElement( rawSplitTree->getParentGlobalIndex(segment) ), stLowerBound );
            }
        }

        // take the intersection and then store it in all corresponding lists.
        T lowerBound = max( jtLowerBound, stLowerBound );
        T upperBound = min( jtUpperBound, stUpperBound );

        for( auto jtSegment : jtContributingSegments ){
            if( !contourTreeSegmentationBounds->count(jtSegment) ){
                (*contourTreeSegmentationBounds)[jtSegment] = unordered_map<int, pair<T,T>>();
            }
            for( auto stSegment : stContributingSegments ){
                (*contourTreeSegmentationBounds)[jtSegment][stSegment] = {lowerBound, upperBound};
                // cout << "writing " << jtSegment << " " << stSegment << " " << lowerBound << " " << upperBound << endl;
            }
        }

        // cout << lowerBound << " " << upperBound << endl;

        T normalLower = sf->getElement(jtContributingSegments[0]);
        T normalUpper = sf->getElement(stContributingSegments[0]);

        if( rawJoinTree->getNode(jtContributingSegments[0])->parent ){
            normalUpper = min(normalUpper, sf->getElement(rawJoinTree->getNode(jtContributingSegments[0])->globalIndex));
        }
        if( rawSplitTree->getNode(stContributingSegments[0])->parent ){
            normalLower = max(normalLower, sf->getElement(rawSplitTree->getNode(stContributingSegments[0])->globalIndex));
        }

        // cout << "(" << normalLower << "," << normalUpper << ") (" << lowerBound << "," << upperBound << ")" << endl;

        // connect the branch in the contour tree
        contourTree->connectNodes( nextLeaf, nextParent );

        // remove the leaf from the join tree and split tree
        // track how the new segmentation merges using the segmentation lists.

        joinTree.deleteNode(nextLeaf, &augmentedJTSegmentation);
        splitTree.deleteNode(nextLeaf, &augmentedSTSegmentation);
        
        // if either parent is now a leaf, add them

        if( nextLeaf == jtParent ){
            exit(0);
        }

        if( jtParent != -1 && (joinTree.getNode(jtParent)->children.size() + splitTree.getNode(jtParent)->children.size() == 1) && queueMembers.find(jtParent) == queueMembers.end() ){
            leafQueue.push_back(jtParent);
            queueMembers.insert(jtParent);
        }
        if( stParent != -1 && (joinTree.getNode(stParent)->children.size() + splitTree.getNode(stParent)->children.size() == 1) && queueMembers.find(stParent) == queueMembers.end() ){
            leafQueue.push_back(stParent);
            queueMembers.insert(stParent);
        }

    }
    // cout << "done" << endl;

    if( verbose ){
        cout << "combined trees" << endl;
    }

    // delete heaps as they will not be needed any more.
    for( auto m : (*minsSorted) ){
        delete m->heap;
    }

    for( auto m : (*maxsSorted) ){
        delete m->heap;
    }

    // clear out these lists
    for( auto m : mins ){
        m->parent = nullptr;
        m->children.clear();
    }

    for( auto m : maxs ){
        m->parent = nullptr;
        m->children.clear();
    }

    MergeTreeInfo<T>* joinTreeInfo = new MergeTreeInfo<T>( rawJoinTree, minsSorted, minSegmentation, jtVTIS, jtUSSC, jtUSCI );
    MergeTreeInfo<T>* splitTreeInfo = new MergeTreeInfo<T>( rawSplitTree, maxsSorted, maxSegmentation, stVTIS, stUSSC, stUSCI );

    return ContourTreeInfo<T>( contourTree, contourTreeNodeIndices, joinTreeInfo, splitTreeInfo, contourTreeSegmentationBounds, unsimplifiedSize );

} // end function