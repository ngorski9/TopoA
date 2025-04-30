#include <iostream>
#include <utility>

#include <fstream>
#include <memory>

#include "scalarfield.hpp"
#include "mergeTreeHelper.hpp"
#include "treeInfoClasses.hpp"

template <typename T> 
void doIt(string file, string outFile, bool verbose){
    ScalarField<T> sf = ScalarField<T>();
    sf.loadFromVTK(file, "Scalars_");
    pair<MergeTreeInfo<T>*,MergeTreeInfo<T>*> trees = getJoinAndSplitTrees(&sf, (T)0, verbose);

    int jtHeight = 0;
    int stHeight = 0;

    for(int which = 0; which < 2; ++which){
        MergeTreeInfo<T>*& tree = (which == 0) ? trees.first : trees.second;
        int& height = (which == 0) ? jtHeight : stHeight;

        int numNodes = tree->vertices->size();

        // identify root:

        std::shared_ptr<CriticalPoint<T>> root = (*(tree->vertices))[0];
        while(root->parent){
            root = root->parent;
        }

        std::vector<pair<std::shared_ptr<CriticalPoint<T>>,int>> pointsQueue;
        pointsQueue.reserve(tree->tree->numNodes);
        pointsQueue.push_back({root,0});

        for(int i = 0; i < pointsQueue.size(); ++i){
            pair<std::shared_ptr<CriticalPoint<T>>,int> nextPair = pointsQueue[i];
            if( nextPair.second > height ){
                height = nextPair.second;
            }
            for ( auto c : nextPair.first->children ){
                pointsQueue.push_back({c, nextPair.second+1});
            }
            cout << i << "/" << numNodes << "                   \r";
        }
        cout << "                                          " << endl;
    }

    ofstream out;
    out.open(outFile, std::ios_base::app);
    out << file << "," << jtHeight << "," << stHeight << endl;

}

int main(int argc, char* argv[]){
    vector<string> args;

    for(int i = 1; i < argc; ++i){
        args.push_back(argv[i]);
    }

    if(argc < 3){
        std::cout << "ERROR: usage is [file] [outputFile] (--verbose) (--singlePrecision) " << std::endl;
        exit(1);
    }

    bool verbose = false;
    bool singlePrecision = false;

    for( int i = 2; i < args.size(); ++i ){
        if( args[i] == "--verbose" ){
            verbose = true;
        } else if( args[i] == "--singlePrecision" ){
            singlePrecision = true;
        } else {
            std::cout << "unrecognoized parameter " << args[i] << std::endl;
            exit(1);
        }
    }

    if( singlePrecision ){
        doIt<float>(args[0], args[1], verbose);
    } else {
        doIt<double>(args[0], args[1], verbose);
    }

}