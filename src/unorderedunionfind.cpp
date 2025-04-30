#include "unorderedunionfind.hpp"
#include <iostream>
void UnorderedUnionFind::add_(int idx){
    if( entries.find(idx) == entries.end() ){
        entries[idx] = idx;
        originalTreeTracker[idx] = idx;
    }
}

const int studying = -1;

// removes all links between the child and parent
void UnorderedUnionFind::resetBetween_(int child, int parent){

    int topOfChain = find_(parent);

    while( child != parent ){

        int nextParent = originalTreeTracker[child];
        if( nextParent == child ){
            break;
        }

        vector<int> resetToChild;
        resetToChild.push_back(child);

        // discover all elements that that descend from the child.
        for( int i = 0; i < resetToChild.size(); ++i ){
            
            if( reverseTracking.count(resetToChild[i]) ){
                resetToChild.reserve( resetToChild.size() + reverseTracking[resetToChild[i]].size() );
                for( auto childOfResetNode : reverseTracking[resetToChild[i]] ){
                    resetToChild.push_back(childOfResetNode);
                }
            }

        }

        for( int i = 0; i < resetToChild.size(); ++i ){
            if( find_( resetToChild[i] ) == topOfChain ){
                entries[resetToChild[i]] = child;
            }
        }


        originalTreeTracker[child] = child;
        vector<int>& parentReverseTracking = reverseTracking[nextParent];
        parentReverseTracking.erase( find( parentReverseTracking.begin(), parentReverseTracking.end(), child ) );

        child = nextParent;

    }
}

// can also lead to buggy behavior. Be careful!
void UnorderedUnionFind::reset_(int idx){

    if( idx == 3511524 ){
        cout << "resetting 3511524" << endl;
    }

    if( idx == studying ){
        cout << "reset studying" << endl;
    }
    int topOfChain = find_(idx);    
    // cout << "d" << endl;
    // go through every node that merges into the current node, and if it still does
    // (i.e., no lower node has been removed from the tree), then link it to the
    // highest possible node before reaching the one that is reset
    for( auto child : reverseTracking[idx] ){
        // cout << "child " << child << endl;
        if( child == idx ){
            continue;
        }
        vector<int> resetToChild;
        // unordered_set<int> elementsOfResetToChild;
        resetToChild.push_back(child);

        // discover all elements that that descend from the child.
        for( int i = 0; i < resetToChild.size(); ++i ){
            
            if( reverseTracking.count(resetToChild[i]) ){
                resetToChild.reserve( resetToChild.size() + reverseTracking[resetToChild[i]].size() );
                for( auto childOfResetNode : reverseTracking[resetToChild[i]] ){
                    // if( elementsOfResetToChild.count(childOfResetNode) == 0 ){
                        resetToChild.push_back(childOfResetNode);
                    //     elementsOfResetToChild.insert(childOfResetNode);
                    // }
                }
            }

        }

        // if there is a path from these elements to the parent of the node that we are removing
        // then cut off that chain at the child.
        for( int i = 0; i < resetToChild.size(); ++i ){
            if( find_( resetToChild[i] ) == topOfChain ){
                entries[resetToChild[i]] = child;
            }
        }

        // for the child itself, say that its original parent in the tree is itself (i.e. not the node above it)
        originalTreeTracker[child] = child;

    }

    // cout << "c" << endl;
    if( originalTreeTracker[idx] != idx ){
        // cout << "i" << endl;
        vector<int>& parentReverseTracking = reverseTracking[originalTreeTracker[idx]];
        // cout << "ii" << endl;
        if( find( parentReverseTracking.begin(), parentReverseTracking.end(), idx ) == parentReverseTracking.end() ){
            cout << idx << endl;
            cout << "it is not present" << endl;
            exit(1);
        }
        parentReverseTracking.erase( find( parentReverseTracking.begin(), parentReverseTracking.end(), idx ) );        
        // cout << "iii" << endl;
        originalTreeTracker[idx] = idx;
    }
    // cout << "b" << endl;

    entries[idx] = idx;
    reverseTracking[idx].clear();
    // cout << "a" << endl;

}

void UnorderedUnionFind::union_(int idx1, int idx2){

    if( find_(idx1) != idx2 && find_(idx1) != find_(idx2) ){
        if( reverseTracking.count(idx2) == 0){
            reverseTracking[idx2] = vector<int>();
            reverseTracking[idx2].reserve(2);
        }
        reverseTracking[idx2].push_back(find_(idx1));

        if( find_(idx1) == studying ){
            cout << "linking studying to " << idx2 << endl;
            cout << "currently, it is linked to " << find_(idx1) << endl;
            cout << "reverse of parent:" << endl;
            for( auto i : reverseTracking[idx2] ){
                cout << "\t" << i << endl;
            }
            cout << "reverse of self:" << endl;
            for( auto i : reverseTracking[originalTreeTracker[idx1]] ){
                cout << "\t" << i << endl;
            }
        }

        if( idx2 == studying ){
            cout << "linking " << find_(idx1) << " to studying" << endl;
        }

        // cout << "connecting " << find_(idx1) << " to " << idx2 << endl;
        originalTreeTracker[find_(idx1)] = idx2;
        entries[find_(idx1)] = idx2;

    }
}

int UnorderedUnionFind::find_(int idx){
    int next_idx = entries[idx];
    int temp = -1;

    while( idx != next_idx ){
        temp = next_idx;
        next_idx = entries[next_idx];
        entries[idx] = next_idx;
        idx = temp;
    }

    return idx;
}