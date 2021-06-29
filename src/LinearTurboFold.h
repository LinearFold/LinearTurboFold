/*
 *LinearTurboFold.h*
 header file for LinearTurboFold.cpp.

 author: Sizhen Li
 edited by: 10/2020
*/

#ifndef LINEAR_TURBOFOLD_H
#define LINEAR_TURBOFOLD_H

#include <string>
#include <vector>

#include "ConfigParser.h"
#include "utils/TProgressDialog.h"
#include "utils/MultiSequence.h"
#include "utils/GuideTree.h"
#include "SeqFold.h"
#include "utils/phmm_aln.h"
#include "utils/phmm.h"
#include "utils/p_alignment.h"
#include "LinearAlignment/src/LinearAlign.h"
#include "ProbabilisticModel.h"
#include "Alignment.h"

enum{NO_ERRORS,
CONSTRUCTOR_ERROR, 
ISEQ_ARGUMENT_OVERFLOW_ERROR,
THREAD_SCHEDULE_ERROR,
UNFINISHED_THREAD_ERROR,
RUN_REFOLDING_THREAD_ERROR,
ALIGNMENT_INFO_COMPUTATION_ERROR,
RNALIB_PROBKNOT_ERROR, 
RNALIB_THRESHOLDING_ERROR, 
RNALIB_MEA_ERROR,
RNALIB_PARTITIONFUNCTION_ERROR,
RNALIB_GETPAIR_ERROR,
RNALIB_WRITECT_ERROR,
RNALIB_GETPAIRPROBABILITY_ERROR,
RNALIB_READSHAPE_ERROR,
RNALIB_SETTEMPERATURE_ERROR,
RNA_LIB_RNA_CONSTRUCTOR_ERROR,
RNALIB_SETDISTANCE_ERROR,
N_TF_ERRORS};

class LinearTurboFold {
    public:
        int err_code;//Internal error code
        int num_seq;
        int n_iterations;
        ProgressHandler *progress;//Used to provide progress information to the interface
        MultiSequence *multiple_sequences;
	    MultiSequence *multiple_alignment;
        string lastErrorDetails;
        vector<t_structure*> sequences;
        vector<SeqFold*> folds;
        vector<char*> saves;
        string aln_save_name;

        BeamAlign* linearAlign;
        double** similarities; // Similarities between the sequences, between 0 and 1.
        vector<vector<unordered_map<int, aln_ret>*>> beam_aln_results;
        vector<unordered_map<int, match_score>> match_score_results;
        match_score** match_score_rets;

        //beam size
        int hmm_beam = 0;
        int cky_beam = 0;
        bool isVerbose = false;

        LinearTurboFold(vector<t_structure*> *fasta_sequences, int ckyBeam, int hmmBeam, bool verbose);
        ~LinearTurboFold();

        void run_iterations(ConfigParser config);
        void SetProgress(ProgressHandler& Progress);
        int GetNumberSequences();

        int allocate_phmm();
        void initialize_multiple_sequences();
        int run_phmm_alignment(int i_iter);
        int run_multiple_alignment();
        double get_match_score(int i_seq1, int i_seq2, int i, int k);
        
        void refoldSequences(int i_iter, bool saveBpps, bool savePfs, vector<string>* pfSaveFiles, vector<string>* bppSaveFiles);
        double get_folding_extrinsic_information(int i_seq1, int j, int i);
        int ProbKnot(const int i_seq, const int n_iterations, const int minhelixlength);
        int WriteCt(const int i_seq, const char fp[]);
};




#endif // LINEAR_TURBOFOLD_H