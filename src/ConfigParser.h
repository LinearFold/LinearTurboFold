/*
 *ConfigParser.h*
 header file for ConfigParser.cpp.

 author: Sizhen Li
 edited by: 10/2020
*/

#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <string>
#include <vector>
// #include "utils/ParseCommandLine.h"
#include "utils/common_utils.h"
#include "utils/structure_object.h"

class ConfigParser {
    public:
        string mode;
        string OutAln;
        string conf_file;
        string AlnFormat;
        double percent;
        double meaGamma;
        double threshold;
        double turboGamma;
        double temperature;
        double turboIterations;
        int windowSize;
        int ColumnNumber;
        int pkIterations;
        int maxStructures;
        int minHelixLength;
        vector<string> sequenceFiles, outputCtFiles, outputDotFiles, outputPfsFiles, shapeFiles, outputBppFiles;
        vector<t_structure*> *fasta_sequences;

        ConfigParser();
        void ParseConfig(int argc, char** argv);

        // linearturbofold
        int hmmbeam;
        int ckybeam;
        bool verbose;
        bool saveBpps;
        bool savePfs;
};

#endif