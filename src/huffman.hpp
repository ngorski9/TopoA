#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <utility>

using namespace std;

class HuffmanNode{
public:
    vector<char> byteCode; // initially, stored backwards because that makes it easier to calculate.
    HuffmanNode* leftChild;
    HuffmanNode* rightChild;
    int weight;
    int symbol;

    HuffmanNode(int symbol_, int weight_){
        weight = weight_;
        symbol = symbol_;
        leftChild = NULL;
        rightChild = NULL;
    }

    HuffmanNode(HuffmanNode* child1, HuffmanNode* child2){
        leftChild = child1;
        rightChild = child2;
        weight = child1->weight + child2->weight;
    }

    ~HuffmanNode(){
        delete leftChild;
        delete rightChild;
    }
};

// greater because make heap works in reverse
class HuffmanNodeCompare{
    public:
        bool operator()(HuffmanNode*& left, HuffmanNode*& right){
            return (left->weight > right->weight);
        }
};

// returns a pair of the root node and a list of all of the leaf nodes.
HuffmanNode* makeHuffmanTree( vector<int>& keys, unordered_map<int,int>& freqs ){

    vector<HuffmanNode*> treeNodesHeap;
    treeNodesHeap.reserve(keys.size());

    for( auto k : keys ){
        HuffmanNode* nextNode = new HuffmanNode( k, freqs[k] );
        treeNodesHeap.push_back(nextNode);
    }

    HuffmanNodeCompare compare;

    make_heap(treeNodesHeap.begin(), treeNodesHeap.end(), compare);

    while( treeNodesHeap.size() > 1 ){

        HuffmanNode* child1 = treeNodesHeap[0];
        pop_heap(treeNodesHeap.begin(), treeNodesHeap.end(), compare);
        treeNodesHeap.pop_back();

        HuffmanNode* child2 = treeNodesHeap[0];
        pop_heap(treeNodesHeap.begin(), treeNodesHeap.end(), compare);
        treeNodesHeap.pop_back();

        HuffmanNode* parent = new HuffmanNode(child1, child2);

        treeNodesHeap.push_back(parent);
        push_heap(treeNodesHeap.begin(), treeNodesHeap.end(), compare);

    }

    return treeNodesHeap[0];

}

void recursiveGetHuffmanCodes( unordered_map<int, vector<char>>& codes, vector<char>& codeSoFar, HuffmanNode* node ){
    if( node->leftChild == NULL ){
        codes[node->symbol] = codeSoFar;
    } else {
        codeSoFar.push_back(0);
        recursiveGetHuffmanCodes( codes, codeSoFar, node->leftChild );
        codeSoFar.pop_back();
        codeSoFar.push_back(1);
        recursiveGetHuffmanCodes( codes, codeSoFar, node->rightChild );
        codeSoFar.pop_back();
    }
}

// first is the table header, second is the list of bytes
pair<vector<int>, vector<char>> huffmanEncode( vector<int> nums ){

    vector<int> keys;
    unordered_map<int, int> freqs;

    for( auto k : nums ){
        if( freqs.count(k) ){
            ++freqs[k];
        } else {
            freqs[k] = 1;
            keys.push_back(k);
        }
    }

    vector<int> tableHeader;
    tableHeader.reserve(2*keys.size());
    for( auto k : keys ){
        tableHeader.push_back(k);
        tableHeader.push_back(freqs[k]);
    }

    HuffmanNode* huffmanTree = makeHuffmanTree(keys, freqs);

    unordered_map<int, vector<char>> codes;
    vector<char> tempVector;
    recursiveGetHuffmanCodes(codes, tempVector, huffmanTree);

    vector<char> encoded;
    char nextByte = 0;
    int nextPosition = 0;

    for( auto n : nums ){
        for( auto bit : codes[n] ){
            if( bit == 1 ){
                nextByte |= (1 << nextPosition);
            }

            ++nextPosition;
            if( nextPosition == 8 ){
                encoded.push_back(nextByte);
                nextByte = 0;
                nextPosition = 0;
            }
        }
    }

    if( nextPosition != 0 ){
        encoded.push_back(nextByte);
    }

    return {tableHeader, encoded};

}

// used for streaming out the huffman numbers.
class HuffmanDecoder{

    HuffmanNode* huffmanTree;
    char* inputCode;
    int nextByte = 0;
    int nextBit = 0;

public:
    HuffmanDecoder( vector<int>& tableHeader, char* inputCode_ ){
        vector<int> keys;
        unordered_map<int,int> freqs;

        int numKeys = tableHeader.size() / 2;
        keys.reserve(numKeys);
        for( int i = 0; i < numKeys; ++i ){
            keys.push_back(tableHeader[2*i]);
            freqs[tableHeader[2*i]] = tableHeader[2*i+1];
        }

        huffmanTree = makeHuffmanTree(keys, freqs);
        inputCode = inputCode_;
    }

    ~HuffmanDecoder(){
        delete huffmanTree;
    }

    int nextSymbol(){
        HuffmanNode* iter = huffmanTree;

        // checking for the left child is the same as checking if it is an interior node.
        while( iter->leftChild ){
            if( inputCode[nextByte] & (1 << nextBit) ){
                iter = iter->rightChild;
            } else {
                iter = iter->leftChild;
            }

            ++nextBit;
            if( nextBit == 8 ){
                nextBit = 0;
                ++nextByte;
            }
        }

        return iter->symbol;
    }

};