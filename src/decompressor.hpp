#include "huffman.hpp"
#include "scalarfield.hpp"
#include "results.hpp"
#include "cubicsplineinterpolation.hpp"
#include "compressorData.hpp"
#include "python_error.hpp"

#include <chrono>
#include <cstdlib>
#include <string>
#include <filesystem>
#include <unistd.h>

using namespace std;

template <typename T>
void decompress(string outputFolder, string inputFilename, string outputFilename, string baseCompressor, bool verbose,  bool saveToVtk, Results* results, string baseCompressorFolder = "../../base_compressors"){
    string systemSuffix;

    if( !verbose ){
        systemSuffix = " > /dev/null 2> /dev/null";
    }

    auto startTime = std::chrono::high_resolution_clock::now();

    int result = system( ("rm " + outputFolder + "/compressed.tar" + systemSuffix).c_str() );

    if( inputFilename != "" ){
        result = system( ("cp " + inputFilename + " " + outputFolder + "/compressed.tar.xz" + systemSuffix).c_str() );
    }

    result = system( ("rm " + outputFolder + "/values.bytes" + systemSuffix).c_str() );
    result = system( ("rm " + outputFolder + "/rawData.cmp" + systemSuffix).c_str() );
    result = system( ("rm " + outputFolder + "/intermediateData.dat" + systemSuffix).c_str() );
    result = system( ("rm " + outputFilename + systemSuffix).c_str() );

    result = system( ("xz -dv " + outputFolder + "/compressed.tar.xz" + systemSuffix).c_str() );

    string cwd = filesystem::current_path();
    result = chdir(outputFolder.c_str());
    result = system( ("tar xvf compressed.tar" + systemSuffix).c_str() );
    result = chdir(cwd.c_str());
    result = system( ("rm " + outputFolder + "/compressed.tar" + systemSuffix).c_str() );

    ifstream in( (outputFolder + "/values.bytes").c_str(), ios::binary );

    int size_x;
    int size_y;
    int size_z;
    T xi;
    int lossless32Size;
    int lossless64Size;
    int tableHeaderSize;
    int huffmanSize;
    double compressorParameter;

    in.read( reinterpret_cast<char*>(&size_x), sizeof(int) );
    in.read( reinterpret_cast<char*>(&size_y), sizeof(int) );
    in.read( reinterpret_cast<char*>(&size_z), sizeof(int) );
    in.read( reinterpret_cast<char*>(&xi), sizeof(T) );
    in.read( reinterpret_cast<char*>(&lossless32Size), sizeof(int) );
    in.read( reinterpret_cast<char*>(&lossless64Size), sizeof(int) );
    in.read( reinterpret_cast<char*>(&tableHeaderSize), sizeof(int) );    
    in.read( reinterpret_cast<char*>(&huffmanSize), sizeof(int) );
    in.read( reinterpret_cast<char*>(&compressorParameter), sizeof(double) );

    int numPoints = size_x * size_y * size_z;

    vector<float> lossless32(lossless32Size);
    vector<double> lossless64(lossless64Size);
    vector<int> tableHeader(tableHeaderSize);
    vector<char> huffmanCodes(huffmanSize);
    vector<int> quantizationCodes;

    in.read( reinterpret_cast<char*>(lossless32.data()), lossless32Size*4 );
    in.read( reinterpret_cast<char*>(lossless64.data()), lossless64Size*8 );
    in.read( reinterpret_cast<char*>(tableHeader.data()), tableHeaderSize*4 );
    in.read( reinterpret_cast<char*>(huffmanCodes.data()), huffmanSize );

    in.close();

    string typeFlag;

    if( baseCompressor == "SZ3" ){

        if( is_same<T, double>::value ){
            typeFlag = "-d";
        } else if( is_same<T, float>::value ){
            typeFlag = "-f";
        }
        result = system( (baseCompressorFolder + "/sz3 " + typeFlag + " -z " + outputFolder + "/rawData.cmp -o " + outputFolder + "/intermediateData.dat -3 " + to_string(size_x) + " " + to_string(size_y) + " " + to_string(size_z) + " -M ABS " + to_string(compressorParameter*xi) + systemSuffix).c_str() );
    } else if( baseCompressor == "ZFP" ){

        if( is_same<T, double>::value ){
            typeFlag = "-d";
        } else if( is_same<T, float>::value ){
            typeFlag = "-f";
        }

        result = system( (baseCompressorFolder + "/zfp " + typeFlag + " -z " + outputFolder + "/rawData.cmp -o " + outputFolder + "/intermediateData.dat -3 " + to_string(size_x) + " " + to_string(size_y) + " " + to_string(size_z) + " -a " + to_string(compressorParameter*xi) + systemSuffix).c_str() );
    } else if( baseCompressor == "TTHRESH" ){

        result = system( (baseCompressorFolder + "/tthresh -c " + outputFolder + "/rawData.cmp -o " + outputFolder + "/intermediateData.dat" + systemSuffix).c_str() );

    } else if( baseCompressor == "CSI" ){
        int ratio = lround(compressorParameter);
        ScalarField<float>* reconstructed32 = csiDecompress<float>(outputFolder + "/rawData.cmp", size_x, size_y, size_z, ratio);
        if( is_same<T,float>::value ){
            reconstructed32->saveToDat(outputFolder + "/intermediateData.dat");
            delete reconstructed32;
        } else {
            ScalarField<double>* reconstructed64 = reconstructed32->upSample();
            reconstructed64->saveToDat(outputFolder + "/intermediateData.dat");
            delete reconstructed32;
            delete reconstructed64;
        }
    } else if( baseCompressor == "Neurcomp" ){

        result = system( ("python3 " + baseCompressorFolder + "/neurcomp/net_decompress.py --volume " + outputFolder + "/dummy.npy --compressed " + outputFolder + "/rawData.cmp --recon " + outputFolder + "/recon" + systemSuffix).c_str() );
        pythonError(result, verbose);

        ScalarField<T> reconVTI;
        reconVTI.loadFromVTK(outputFolder + "/recon.vti", "sf");
        reconVTI.saveToDat(outputFolder + "/intermediateData.dat");
        result = system( ("rm " + outputFolder + "/recon.vti").c_str() );

    }

    ScalarField<T> sf;
    sf.loadFromDat(outputFolder + "/intermediateData.dat", size_x, size_y, size_z);

    HuffmanDecoder* hd = NULL;
    hd = new HuffmanDecoder( tableHeader, huffmanCodes.data() );

    int lossless32Symbol = lround( pow( 2, LOSSLESS_PRECISION_32 + 5 ) ) + 1;
    int lossless64Symbol = lround( pow( 2, LOSSLESS_PRECISION_64 + 5) ) + 1;
    T codeMultiplier = xi * pow( 2, -MAX_PRECISION+1 );

    int nextLossless32 = 0;
    int nextLossless64 = 0;

    int lastSymbol = 0;
    for( int i = 0; i < numPoints; ++i ){
        int nextSymbol;
        nextSymbol = hd->nextSymbol();

        if( nextSymbol == lossless32Symbol ){
            sf.setElement(i, lossless32[nextLossless32]);
            ++nextLossless32;
        } else if( nextSymbol == lossless64Symbol ){
            sf.setElement(i, lossless64[nextLossless64]);
            ++nextLossless64;
        } else {
            sf.setElement(i, sf.getElement(i) + nextSymbol * codeMultiplier );
        }
    }

    if( saveToVtk ){
        sf.saveToVTK(outputFilename);

    } else {
        sf.saveToDat(outputFilename);
    }

    result = system(("rm " + outputFolder + "/" + inputFilename + ".tar" + systemSuffix).c_str());
    result = system(("rm " + outputFolder + "/intermediateData.dat" + systemSuffix).c_str());
    result = system(("rm " + outputFolder + "/rawData.cmp " + systemSuffix).c_str());
    result = system(("rm " + outputFolder + "/values.bytes" + systemSuffix).c_str());

    auto endTime = std::chrono::high_resolution_clock::now();
    double dif = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime).count();
    cout << "decompression time: " << dif << endl;

    if( outputFilename != outputFolder + "/output.format" ){
        cout << "decompressed file: " << (outputFilename) << endl;
    }

    delete hd;

    if( results ){
        results->decompressionTime = dif;
    }
}