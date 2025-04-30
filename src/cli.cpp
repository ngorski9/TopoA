#include "compressor.hpp"
#include "compressorIterativeCT.hpp"
#include "compressorIterativeMT.hpp"
#include "decompressor.hpp"
#include "evaluation.hpp"
#include "results.hpp"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <sys/stat.h>

using namespace std;

void printHelpMessage();

inline bool file_exists (string fname) {
  struct stat buffer;   
  return stat(fname.c_str(), &buffer) == 0;
}

void printUsage(){
    cout << "Run TopoA -h or TopoA -help for more information." << endl;
    exit(1);
}

void printNotEnoughInfo(int i, int stop, string flag){
    if( i == stop ){
        cout << "TopoA: Please specify a value after option '" + flag + "'." << endl;
        printUsage();
    }
}

int main(int argc, char* argv[]){
    vector<string> args;

    if(argc == 1){
        printHelpMessage();
    }

    for(int i = 1; i < argc; ++i){
        args.push_back(argv[i]);
    }

    string filename;
    bool specifiedFilename = false;

    string array;
    bool specifiedArray = false;

    double epsilon;
    bool specifiedEpsilon = false;

    double xi;
    bool specifiedXi = false;

    string outputFolder = ".";
    bool specifiedOutputFolder = false;

    string baseCompressor;
    bool specifiedBaseCompressor = false;

    double compressorParameter;
    bool specifiedCompressorParameter = false;

    string csv;
    bool specifiedCsv = false;

    bool verbose = false;
    bool doublePrecision = true;

    string compressedName = "";
    string decompressedName = "";
    bool saveToVtk = false; // otherwise it saves to a .raw

    bool iterative = false;
    bool iterativeMergeTree = false;

    bool logQuantization = true;
    int initialPrecision = 0;

    bool experiment = false;
    bool checkCorrectness = true;

    string baseCompressorFolder = ".";

    int i = 0;
    while( i < argc - 1){
        if( args[i] == "-i" ){
            printNotEnoughInfo(i, argc-2, "-i");
            filename = args[i+1];
            specifiedFilename = true;
        } else if( args[i] == "-eps" ){
            printNotEnoughInfo(i, argc-2, "-eps");            
            try{
                epsilon = stod(args[i+1]);
            } catch (std::invalid_argument) {
                cout << "TopoA: Persistence threshold epsilon must be a floating point value." << endl;
                printUsage();
            }
            specifiedEpsilon = true;
        } else if( args[i] == "-xi" ){
            printNotEnoughInfo(i, argc-2, "-xi");
            try{
                xi = stod(args[i+1]);
            } catch (std::invalid_argument) {
                cout << "TopoA: Error bound xi must be a floating point value." << endl;
                printUsage();
            }
            specifiedXi = true;
        } else if( args[i] == "-of" ){
            printNotEnoughInfo(i, argc-2, "-of");            
            outputFolder = args[i+1];
            specifiedOutputFolder = true;
        } else if( args[i] == "-bc" ){
            printNotEnoughInfo(i, argc-2, "-bc");            
            baseCompressor = args[i+1];
            specifiedBaseCompressor = true;
        } else if( args[i] == "-bcFolder" ){
            printNotEnoughInfo(i, argc-2, "-bcFolder");
            baseCompressorFolder = args[i+1];
        } else if( args[i] == "-bcParameter" ){
            printNotEnoughInfo(i, argc-2, "-bcParameter");
            compressorParameter = stod(args[i+1]);
            specifiedCompressorParameter = true;
        } else if( args[i] == "-csv" ){
            printNotEnoughInfo(i, argc-2, "-csv");            
            csv = args[i+1];
            specifiedCsv = true;
        } else if( args[i] == "-a" ){
            printNotEnoughInfo(i, argc-2, "-a");
            array = args[i+1];
            specifiedArray = true;
        } else if( args[i] == "-f" ){
            doublePrecision = false;
            --i;
        } else if( args[i] == "-verbose" ){
            verbose = true;
            --i;
        } else if( args[i] == "-c" ){
            printNotEnoughInfo(i, argc-2, "-c");
            compressedName = args[i+1];
        } else if( args[i] == "-o" ){
            printNotEnoughInfo(i, argc-2, "-o");            
            decompressedName = args[i+1];
        } else if( args[i] == "-vtk" ){ 
            saveToVtk = true;
            --i;
        } else if( args[i] == "-iterative" ){
            iterative = true;
            --i;
        } else if( args[i] == "-iterativeMergeTree" ){
            iterativeMergeTree = true;
            --i;
        } else if( args[i] == "-linearQuantization" ){
            logQuantization = false;
            initialPrecision = 1;
            --i;
        } else if( args[i] == "-experiment" ){
            experiment = true;
            checkCorrectness = true;
            --i;
        } else if( args[i] == "-skipCorrectness" ){
            checkCorrectness = false;
            --i;
        } else if( args[i] == "-help" || args[i] == "-h"){
            printHelpMessage();
        }else {
            cout << "Unrecognized argument " << args[i] << endl;
            printUsage();
        }

        ++i;
        ++i;
    }

    // if parameter is not specified, then use the value that we have found to be approximately the best (so far)
    if( !specifiedCompressorParameter ){
        if( baseCompressor == "SZ3" ){
            compressorParameter = 0.25;
        } else if( baseCompressor == "TTHRESH" ){
            compressorParameter = 0.05;
        } else if( baseCompressor == "ZFP" ){
            compressorParameter = 5;
        } else if( baseCompressor == "Neurcomp" ){
            compressorParameter = 9;
        } else if( baseCompressor == "CSI" ){
            compressorParameter = 7;
        }
    }

    bool bad = false;

    // check correctness of options
    if(!specifiedFilename && compressedName == ""){
        cout << "TopoA: You must specify an input filename to compress (with -i) or a compressed file to decompress (with -c)" << endl;
        bad = true;
    }

    if(specifiedFilename && !file_exists(filename)){
        cout << "TopoA: Cannot stat input filename '" + filename + "'." << endl;
        bad = true;
    }

    if(specifiedFilename && compressedName == "" && decompressedName == "" && !experiment){
        cout << "TopoA: When specifying a file to compress, you must either specify at least one of: a compressed filename (with -c), a decompressed filename (with -o), or the experiment setting (with -experiment)" << endl;
        bad = true;
    }

    if(specifiedFilename && !specifiedArray){
        cout << "TopoA: When compressing a file, you must specify the name of the array in the VTK file that you are compressing" << endl;
        bad = true;
    }

    if(!specifiedFilename && compressedName != "" && decompressedName == ""){
        cout << "TopoA: When specifying a file to decompress, you must provide a name for the decompressed file (with -o)." << endl;
        bad = true;
    }

    if(!specifiedFilename && experiment){
        cout << "TopoA: Cannot use the experiment setting without specifying a file to compress." << endl;
        bad = true;
    }

    if(specifiedFilename && (!specifiedEpsilon || epsilon <= 0.0)){
        cout << "TopoA: When compressing, you must specify a persistence threshold epsilon greater than zero." << endl;
        bad = true;
    }

    if(specifiedFilename && (!specifiedXi || xi <= 0.0)){
        cout << "TopoA: When compressing, you must specify an error bound xi greater than zero." << endl;
        bad = true;
    }

    if( compressorParameter <= 0.0 ){
        cout << "TopoA: Compressor parameter must be greater than zero." << endl;
        bad = true;
    }

    if(!file_exists(outputFolder)){
        cout << "TopoA: Cannot stat output folder '" + outputFolder + "'." << endl;
        bad = true;
    }

    if( baseCompressor != "ZFP" && baseCompressor != "SZ3" && baseCompressor != "CSI" && baseCompressor != "TTHRESH" && baseCompressor != "Neurcomp" ){
        cout << "TopoA: Invalid base compressor '" + baseCompressor + "'." << endl;
        cout << "\tPlease specify one of: ZFP, SZ3, TTHRESH, CSI, Neurcomp" << endl;
        bad = true;
    }

    if(!file_exists(baseCompressorFolder)){
        cout << "TopoA: Cannot stat base compressor folder '" + baseCompressorFolder + "'." << endl;
        bad = true;
    }

    if( baseCompressor == "ZFP" && !file_exists( baseCompressorFolder + "/zfp" ) ){
        cout << "TopoA: Cannot stat ZFP binary at '" + baseCompressorFolder + "/zfp'." << endl;
        bad = true;
    }

    if( baseCompressor == "SZ3" && !file_exists( baseCompressorFolder + "/sz3" ) ){
        cout << "TopoA: Cannot stat SZ3 binary at '" + baseCompressorFolder + "/sz3'." << endl;
        bad = true;
    }

    if( baseCompressor == "TTHRESH" && !file_exists( baseCompressorFolder + "/tthresh" ) ){
        cout << "TopoA: Cannot stat TTHRESH binary at '" + baseCompressorFolder + "/tthresh'." << endl;
        bad = true;
    }

    if( baseCompressor == "CSI" ){
        long param = lround(compressorParameter);
        if( param != compressorParameter  || param < 2 || param > 16){
            cout << "TopoA: When using CSI, parameter must be an integer between 2 and 16 (inclusive)." << endl;
            bad = true;
        }
    }

    if( baseCompressor == "Neurcomp" ){
        if( !file_exists( baseCompressorFolder + "/neurcomp" ) ){
            cout << "TopoA: Cannot stat Neurcomp repo at '" + baseCompressorFolder + "/neurcomp'." << endl;
            bad = true;
        }

        if( !file_exists( outputFolder + "/thenet.pth" ) ){
            cout << "TopoA: Cannot stat '" + outputFolder + "/thenet.pth' required for Neurcomp." << endl;
            bad = true;
        }

        if( !file_exists( outputFolder + "/thenet.json" ) ){
            cout << "TopoA: Cannot stat '" + outputFolder + "/thenet.json' required for Neurcomp." << endl;
            bad = true;
        }

        if( !file_exists( outputFolder + "/dummy.npy" ) ){
            cout << "TopoA: Cannot stat '" + outputFolder + "/dummy.npy' required for Neurcomp." << endl;
            bad = true;
        }

        long param = lround(compressorParameter);
        if( param != compressorParameter ){
            cout << "TopoA: When using Neurcomp, parameter must be an integer." << endl;
            bad = true;
        }
    }

    if( iterative && iterativeMergeTree ){
        cout << "TopoA: Please only specify one of -iterative or -iterativeMergeTree." << endl;
        bad = true;
    }

    if( !checkCorrectness && !experiment ){
        cout << "TopoA: Only use -skipCorrectness with -experiment" << endl;
        bad = true;
    }

    if( specifiedCsv && !experiment ){
        cout << "TopoA: Only use -csv with -experiment" << endl;
        bad = true;
    }

    if( bad ){
        printUsage();
    }

    string systemSuffix = "";
    if( !verbose ){
        systemSuffix = " > /dev/null 2> /dev/null";
    }

    if( doublePrecision ){
        Results* results = NULL;
        
        if( specifiedFilename ){

            if( iterative ){
                results = compressIterative<double>(filename, array, epsilon, xi, compressedName, outputFolder, baseCompressor, compressorParameter, verbose, logQuantization, initialPrecision, baseCompressorFolder);
            } else if( iterativeMergeTree ){
                results = compressIterativeMT<double>(filename, array, epsilon, xi, compressedName, outputFolder, baseCompressor, compressorParameter, verbose, logQuantization, initialPrecision, baseCompressorFolder);
            } else{
                results = compress<double>(filename, array, epsilon, xi, compressedName, outputFolder, baseCompressor, compressorParameter, verbose, logQuantization, initialPrecision, baseCompressorFolder);
            }
        }

        if( experiment || decompressedName != "" ){

            string decompressedNameOutput;
            if( decompressedName == "" ){
                decompressedNameOutput = outputFolder + "/output.format";
            } else {
                decompressedNameOutput = decompressedName;
            }

            string compressedNameInput;
            if( specifiedFilename ){
                compressedNameInput = "";
            } else {
                compressedNameInput = compressedName;
            }

            decompress<double>(outputFolder, compressedNameInput, decompressedNameOutput, baseCompressor, verbose, saveToVtk, results, baseCompressorFolder);

            if( experiment ){
                if( specifiedCsv ){
                    writeEvaluationToCSV<double>( csv, baseCompressor, compressorParameter, filename, array, decompressedNameOutput, saveToVtk, xi, epsilon, results, iterative, iterativeMergeTree, logQuantization, initialPrecision, checkCorrectness );
                } else {
                    printEvaluationToConsole<double>( filename, array, decompressedNameOutput, saveToVtk, xi, epsilon, results, checkCorrectness );
                }   
            }

            if( decompressedName == "" ){
                int result = system(("rm " + outputFolder + "/output.format " + systemSuffix).c_str());
            }

        }

        int result = system(("rm " + outputFolder + "/compressed.tar.xz " + systemSuffix).c_str());

    } else {
        Results* results = NULL;
        
        if( specifiedFilename ){

            if( iterative ){
                results = compressIterative<float>(filename, array, epsilon, xi, compressedName, outputFolder, baseCompressor, compressorParameter, verbose, logQuantization, initialPrecision, baseCompressorFolder);
            } else if( iterativeMergeTree ){
                results = compressIterativeMT<float>(filename, array, epsilon, xi, compressedName, outputFolder, baseCompressor, compressorParameter, verbose, logQuantization, initialPrecision, baseCompressorFolder);
            } else{
                results = compress<float>(filename, array, epsilon, xi, compressedName, outputFolder, baseCompressor, compressorParameter, verbose, logQuantization, initialPrecision, baseCompressorFolder);
            }
        }

        if( experiment || decompressedName != "" ){

            string decompressedNameOutput;
            if( decompressedName == "" ){
                decompressedNameOutput = outputFolder + "/output.format";
            } else {
                decompressedNameOutput = decompressedName;
            }

            string compressedNameInput;
            if( specifiedFilename ){
                compressedNameInput = "";
            } else {
                compressedNameInput = compressedName;
            }

            decompress<float>(outputFolder, compressedNameInput, decompressedNameOutput, baseCompressor, verbose, saveToVtk, results, baseCompressorFolder);

            if( experiment ){
                if( specifiedCsv ){
                    writeEvaluationToCSV<float>( csv, baseCompressor, compressorParameter, filename, array, decompressedNameOutput, saveToVtk, xi, epsilon, results, iterative, iterativeMergeTree, logQuantization, initialPrecision, checkCorrectness );
                } else {
                    printEvaluationToConsole<float>( filename, array, decompressedNameOutput, saveToVtk, xi, epsilon, results, checkCorrectness );
                }   
            }

            if( decompressedName == "" ){
                int result = system(("rm " + outputFolder + "/output.format " + systemSuffix).c_str());
            }

        }

        int result = system(("rm " + outputFolder + "/compressed.tar.xz " + systemSuffix).c_str());

    }

}

void printHelpMessage(){
    cout << "Usage: TopoA <options>" << endl;
    cout << "Please see the README for more details." << endl;

    cout << endl;

    cout << "-bc : Base compressor. Choose one of: ZFP, SZ3, TTHRESH, CSI, or Neurcomp." << endl;

    cout << endl;

    cout << "Specifying Files" << endl;
    cout << "\t-i <path> : VTK file for compression." << endl;
    cout << "\t-c <path> : Compressed output." << endl;
    cout << "\t-o <path> : Decompressed output." << endl;
    cout << "\t-a <array name> : " << endl;
    cout << "\t\tName of the array to compress. Not required if you are" << endl;
    cout << "\t\tonly decompressing." << endl;
    cout << "\t-f :" << endl;
    cout << "\t\tSet this flag if the array data uses single precision." << endl;
    cout << "\t\tIf this flag is not set, double precision will be assumed." << endl;
    cout << "\t-vtk :" << endl;
    cout << "\t\tIf this flag is set, save decompressed file in VTK format." << endl;
    cout << "\t\tIf it is not set, then decompress in raw binary format." << endl;

    cout << endl;

    cout << "Compression Parameters (not required for decompression only)." << endl;
    cout << "\t-eps : Persistence threshold. Percent of function range in (0,1]." << endl;
    cout << "\t-xi : Pointwise error bound. Percent of function range in (0,1]." << endl;

    cout << endl;

    cout << "Options for Running Experiments." << endl;
    cout << "\t-experiment :" << endl;
    cout << "\t\tRun compression and decompression, and evaluate the" << endl;
    cout << "\t\tresults. If used, using -c and -o is optional. If neither -c nor" << endl;
    cout << "\t\t-o are set, then TopoA will only work with temporary files." << endl;
    cout << "\t-skipCorrectness :" << endl;
    cout << "\t\tUse with -experiment. If used, will not check whether " << endl;
    cout << "\t\tthe decompressed contour tree is correct (saves time)." << endl;
    cout << "\t-csv <path> :" << endl;
    cout << "\t\tUse with -experiment. Specify CSV file to save results. If file" << endl;
    cout << "\t\talready exists, the results will be appended to the end of the." << endl;
    cout << "\t\tfile." << endl;

    cout << endl;

    cout << "Other Options:" << endl;
    cout << "\t-bcFolder : Specify the folder that contains the base compressor." << endl;
    cout << "\t-of : Specify the output folder where temporary files will be saved." << endl;
    cout << "\t-verbose : Display debug prints and outputs of system calls." << endl;
    cout << "\t-iterative : Compress using an iterative strategy (like TopoSZ)" << endl;
    cout << "\t-iterativeMergeTree :" << endl;
    cout << "\t\tCompress using an iterative strategy, but iterate on the" << endl;
    cout << "\t\tindividual merge trees, rather than the contour tree." << endl;
    cout << "\t-linearQuantization :" << endl;
    cout << "\t\tUse linear-scaling quantization, instead of logarithmic-" << endl;
    cout << "\t\tscaling quantization." << endl;
    cout << "\t-bcParameter." << endl;
    cout << "\t\tSpecify a parameter used for the base compressor (see paper," << endl;
    cout << "\t\tsection 5.2). Specifies the value k for SZ3, ZFP, and TTHRESH." << endl;
    cout << "\t\tSpecifies s for CSI. Specifies 'cluster bits' for Neurcomp." << endl;
    cout << "\t\tDo not specify this if you are only decompressing." << endl;    
    cout << "\t-help, -h : Display this message" << endl;

    exit(0);
}