#pragma once

#include <vector>
#include <cmath>

#include "criticalpointstorage.hpp"
#include "scalarfield.hpp"

template <typename T> class CriticalPointTree{
    public:
        ScalarField<T>* sf;
        // children are referenced as global indices in the scalar field.
        unordered_map<int, shared_ptr<CriticalPoint<T>>> nodes;
        int numNodes = 0;
    public:
        int getNumNodes(){
            return numNodes;
        }

        shared_ptr<CriticalPoint<T>> getNode(int idx){
            return nodes[idx];
        }

        int getParentGlobalIndex( int idx ){
            if( nodes[idx]->parent != nullptr ){
                return nodes[idx]->parent->globalIndex;
            } else {
                return -1;
            }
        }

        bool isLeaf( int idx ){
            return hasNode(idx) && (nodes[idx]->children.size() == 0);
        }

        CriticalPointTree(ScalarField<T>* sf_){
            sf = sf_;
        }

        // aboveIdx is the node that this should be that node's parent.
        void newNode(int criticalType, int criticalIndex, int x, int y, int z, int aboveIdx = -1 ){
            ++numNodes;
            shared_ptr<CriticalPoint<T>> nextNode = shared_ptr<CriticalPoint<T>>(new CriticalPoint<T>(criticalType, criticalIndex, x, y, z, sf));
            nodes[nextNode->globalIndex] = nextNode;

            if(aboveIdx != -1){
                shared_ptr<CriticalPoint<T>> newChild = nodes[aboveIdx];

                nextNode->parent = newChild->parent;
                newChild->parent = nextNode;
                nextNode->children.push_back(newChild);
                if( nextNode->parent != NULL ){
                    vector<shared_ptr<CriticalPoint<T>>>& parentChildren = nextNode->parent->children;
                    parentChildren[ find( parentChildren.begin(), parentChildren.end(), newChild ) - parentChildren.begin() ] = nextNode;
                }
            }
        }

        // this assumes that the child is not currently connected to anything.
        // otherwise it will lead to buggy behavior.
        void connectNodes(int childIdx, int parentIdx){
            shared_ptr<CriticalPoint<T>> child = nodes[childIdx];
            shared_ptr<CriticalPoint<T>> parent = nodes[parentIdx];
            child->parent = parent;
            parent->children.push_back(child);
        }

        // removes a single leaf node from the tree and connects children to its parent.
        // segmentation 
        void deleteNode(int nodeIdx, unordered_map<int,vector<int>>* segmentation = NULL ){
            --numNodes;
            shared_ptr<CriticalPoint<T>> node = nodes[nodeIdx];
            if( node->parent != NULL ){
                vector<shared_ptr<CriticalPoint<T>>>& parentChildren = node->parent->children;
                parentChildren.erase( find(parentChildren.begin(), parentChildren.end(), node) );

                if( segmentation ){
                    for( auto child : node->children ){
                        parentChildren.push_back(child);
                        child->parent = node->parent;
                        
                        vector<int>& childSegmentation = (*segmentation)[child->globalIndex];
                        for( auto segment : (*segmentation)[nodeIdx] ){
                            if( find(childSegmentation.begin(), childSegmentation.end(), segment) == childSegmentation.end() ){
                                childSegmentation.push_back(segment);
                            }
                        }
                    }

                } else {
                    for( auto child : node->children ){
                        parentChildren.push_back(child);
                        child->parent = node->parent;
                    }
                }

            } else {
                if( segmentation ){
                    for( auto child : node->children ){
                        child->parent = NULL;
                        vector<int>& childSegmentation = (*segmentation)[child->globalIndex];
                        for( auto segment : (*segmentation)[nodeIdx] ){
                            if( find(childSegmentation.begin(), childSegmentation.end(), segment) == childSegmentation.end() ){
                                childSegmentation.push_back(segment);
                            }
                        }        
                    }
                } else {
                    for( auto child : node->children ){
                        child->parent = NULL;
                    }
                }
            }
            
            // delete node;
            node->parent = nullptr;
            node->children.clear();
            nodes.erase(nodes.find(nodeIdx));
        }

        bool hasNode(int idx){
            return (nodes.find(idx) != nodes.end());
        }

        T getEdgeLength(int childIdx){
            return abs(sf->getElement(childIdx) - sf->getElement(nodes[childIdx]->parent->globalIndex));
        }
};