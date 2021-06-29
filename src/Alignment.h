#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <iomanip>


#include "utils/SafeVector.h"
#include "utils/MultiSequence.h"
#include "ProbabilisticModel.h"
#include "utils/GuideTree.h"
#include "utils/random.h"
#include "utils/phmm_aln.h"

#include "LinearAlignment/src/LinearAlign.h"

#define numConsistencyReps 2     // Number of probabilistic consistency transformation.
#define numIterativeRefinementReps 100   // Number of iterative refinement of the multiple alignment.


//! This function computes the consistency transformation for sequence Z by taking two posterior probabilities matrices of alignments between X-Z and Z-Y.
//! For the case that sequence Z's index is larger than that of X.
//! The transformed matrix is added to \param posterior.
void LinearConsistencyTransform(int lengthX, unordered_map<int, aln_ret>* &xz_aln_ret, unordered_map<int, multi_aln_ret>* &zy_aln_ret, unordered_map<int, multi_aln_ret>* &new_xy_ret);

//! This function takes multiple sequences and posterior probability matices to perform three-way probabilistic consistency transformation.
//! Returns new re-estimated alignment score matrices. 
//! The formula is: P'(x[i]-y[j])=(1/|S|)*sum_z_in_S{ sum_k{ P(x[i]-z[k]) * P(z[k]-y[j]) } }
vector<vector<unordered_map<int, multi_aln_ret>*>> LinearMultiConsistencyTransform(MultiSequence *sequences, vector<vector<unordered_map<int, multi_aln_ret>*>> &beam_aln_rets);

//! This function takes two multiple sequence alignments as input.
//! Returns the alignment of the two MultiSequence objects.
MultiSequence *LinearAlignAlignments (MultiSequence *align1, MultiSequence *align2,
                                const vector<vector<unordered_map<int, multi_aln_ret>*>> &multi_aln_results,
                                const ProbabilisticModel &model, int hmmBeam);

//! This function takes guide tree (computed by distance) as input.
//! Returns the aligned sequences corresponding to a node or leaf of a guide tree.
MultiSequence *LinearProcessTree (const TreeNode *tree, MultiSequence *sequences,
                            const vector<vector<unordered_map<int, multi_aln_ret>*>> &multi_aln_results,
                            const ProbabilisticModel &model, int hmmBeam);

//! This function computes the final alignment by calling ProcessTree() and performing iterative refinement.
MultiSequence *LinearComputeFinalAlignment (const TreeNode *tree, MultiSequence *sequences,
                                      const vector<vector<unordered_map<int, multi_aln_ret>*>> &multi_aln_results,
                                      const ProbabilisticModel &model, int hmmBeam);

//! This function performs randomized partitioning iterative refinement. 
//! Taking posterior probability matrices, parameters of probabilistic model, and multiple sequence alignments.
//! Returns a new multiple sequence alignment.
void LinearDoIterativeRefinement (const vector<vector<unordered_map<int, multi_aln_ret>*>> &multi_aln_results,
                            const ProbabilisticModel &model, MultiSequence* &alignment, int i, int hmmBeam);