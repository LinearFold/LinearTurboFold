#ifndef RNA_FOLD_H
#define RNA_FOLD_H

#include "utils/structure.h"
#include "LinearPartition/src/LinearPartition.h"
#include "LinearAlignment/src/LinearAlign.h"

typedef int RNAInputType; // Just to help new programmers find the right constants (defined below)
#define SEQUENCE_STRING 0
#define	FILE_SEQ 2

class LinearTurboFold;
class SeqFold {
    public:
        //Integer to keep track of error codes.
		//These errors result on file i/o problems during construction.
		//The errors can be accessed using GetErrorCode().
        int ErrorCode;

        // Holds error messages that occur when reading files, so additional information can be shown to the users (e.g. in GUI programs where console output is not seen).
		string lastErrorDetails;

        //The primitive class for storing sequence and structure data.
		//Private inheritance to force the user to use the interface provided by RNA.
		//call GetStructure() if you need to access the data
		structure *ct;

        //The following bool is used to indicate whether the partion function arrays have been allocated and therefore need to be deleted.
		bool partitionfunctionallocated;

        BeamCKYParser *linearpf;
        LinearTurboFold *outer;

        unordered_map<int, State> *pfscore;
        unordered_map<int, State> *prepfscore;
		unordered_map<int, ExtValue> *extrisic_info;
		unordered_map<int, ExtValue> *pre_ext_info;
		double previterbi;
		char* saveFile;
        double viterbi_score;

        SeqFold(int ckyBeam, const char filepathOrSequence[], bool verbose);
        ~SeqFold();
        void SetSequenceLabel(const string& label);
        int FileReader(const char filename[], const RNAInputType type);
        structure *GetStructure();
        int GetErrorCode() const;
        // string GetFullErrorMessage() const;
        int LinearPartition(int i_iter, int i_seq, string pfsSaveFile, string bppSaveFile);
        void Savepartitionfunction(int i_iter);
        void ComputeMatchScore(match_score** &match_score_ret, int i_seq);
        double GetBasePairProb(int i, int j);

        int ProbKnot(int iterations=1, int MinHelixLength = 1, double threshold = 0);
        int WriteCt(const char filename[], bool append=false, 
			CTCommentProvider &commentProvider=CTComments::Energy) const;
};

#endif