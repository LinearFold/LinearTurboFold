#ifndef _TURBOFOLD_PROBABILISTICMODEL_
#define _TURBOFOLD_PROBABILISTICMODEL_

#include <list>
#include <cmath>
#include <cstdio>
#include <unordered_map>
#include "utils/SafeVector.h"
#include "utils/MultiSequence.h"
#include "LinearAlignment/src/LinearAlign.h"

#ifdef TURBOHOMOLOGY
	#define ChooseBestOfThree ChooseBestOfThreeturbohomology
	#define ProbabilisticModel ProbabilisticModelturbohomology
#endif


using namespace std;

//! ProbabilisticModel Class
/*!
    The ProbabilisticModel Class stores the parameters of a probabilistic model.

*/

#define NumMatrixTypes 3   // One match state and two insert states.

//! Store the largest of three values x1, x2, and x3 in *x.  
//! If xi is the largest value, then store bi in *b.
void ChooseBestOfThree(float x1, float x2, float x3, char b1, char b2, char b3, float *x, char *b);

class ProbabilisticModel{

public:
    //! Computes an alignment based on given posterior matrix.
    //! This is done by finding the maximum summing path (or maximum weight trace) through the posterior matrix.
    //! The final alignment is returned as a pair consisting of:
    //! (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and denote insertions in one of the two sequences and B's denote that both sequences are present (i.e. matches).
    //! (2) a float indicating the sum achieved.
    SafeVector<float> *LinearBuildPosterior(MultiSequence *align1, MultiSequence *align2,
                      const vector<vector<unordered_map<int, multi_aln_ret>*>> &mul_aln_rets, float cutoff = 0.0f) const;

    unordered_map<int, multi_aln_ret> *LinearMultiAlnResults(MultiSequence *align1, MultiSequence *align2,
                   const vector<vector<unordered_map<int, multi_aln_ret>*>> &mul_aln_rets, float cutoff = 0.0f) const;

    pair<SafeVector<char> *, float> LinearComputeAlignment(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, multi_aln_ret>* &mul_aln_rets) const;
};

#endif