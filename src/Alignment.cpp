#include "Alignment.h"

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
#include <unordered_map>


void LinearConsistencyTransform(int lengthX, unordered_map<int, multi_aln_ret>* &xz_aln_ret, unordered_map<int, multi_aln_ret>* &zy_aln_ret, unordered_map<int, multi_aln_ret>* &new_xy_ret){
    // int lengthX = xz_aln_ret->size();
    // cout << "lengthX: " << lengthX << endl;
    // int lengthZ = zy_aln_ret->size();

    for(int i = 1; i <= lengthX; i++){
        // cout << i << " ";  
        for(auto &candXZ : xz_aln_ret[i]){
            int j = candXZ.first;

            for(auto &candZY : zy_aln_ret[j]){
                int k = candZY.first;
                new_xy_ret[i][k].value += candXZ.second.value * candZY.second.value;
                // cout << candXZ.second.aln_prob * candZY.second.aln_prob << " ";
            }
        }
        // cout << endl;
    }
}

// linearTurboFold
vector<vector<unordered_map<int, multi_aln_ret>*>> LinearMultiConsistencyTransform(MultiSequence *sequences, vector<vector<unordered_map<int, multi_aln_ret>*>> &beam_aln_rets){
    const int numSeqs = sequences->GetNumSequences();
    vector<vector<unordered_map<int, multi_aln_ret>*>> new_aln_results;
    new_aln_results.resize(numSeqs);
    for (int i = 0; i < numSeqs; i++){
        new_aln_results[i].resize(numSeqs);
        for (int j = 0; j < numSeqs; j++){
            if (i == j) continue;
            new_aln_results[i][j] = new unordered_map<int, multi_aln_ret>[sequences->GetSequence(i)->GetLength() + 1];
        }
    }

    // For every pair of sequences
    for (int i = 0; i < numSeqs; i++){
        new_aln_results[i].resize(numSeqs);

        for (int j = i+1; j < numSeqs; j++){
            Sequence *seq1 = sequences->GetSequence(i);
            Sequence *seq2 = sequences->GetSequence(j);

            const int seq1Length = seq1->GetLength();
            const int seq2Length = seq2->GetLength();

            // allocate space for temporary results
            unordered_map<int, multi_aln_ret>* tmp_aln_ret = new unordered_map<int, multi_aln_ret>[seq1Length + 1];

            // Get the original alignment result
            unordered_map<int, multi_aln_ret>* &aln_ret = beam_aln_rets[i][j];

            // Contribution from the summation where z = x and z = y
            // cout << "seq1length: " << seq1Length << endl;
            for (int k = 0; k <= seq1Length; k++){
                for(auto &item : aln_ret[k]){
                    int l = item.first;
                    tmp_aln_ret[k][l].value = 2 * aln_ret[k][l].value;
                    // cout << i << " " << j << " " << k << " " << l << " " << new_aln_ret[k][l].value << endl;
                }
            }

            // Contribution from all other sequences
            for (int k = 0; k < numSeqs; k++) {
                if (k == i || k == j) continue;
                // cout << "seqs: " << i << " " << j << " " << k << endl;
                LinearConsistencyTransform(seq1Length, beam_aln_rets[i][k], beam_aln_rets[k][j], tmp_aln_ret);
            }

            // Renormalization
            for (int k = 0; k <= seq1Length; k++){
                for(auto &item : tmp_aln_ret[k]){
                    int l = item.first;
                    tmp_aln_ret[k][l].value /= numSeqs;
                }
            }

            // Mask out positions not originally in the posterior matrix
            for (int k = 0; k <= seq1Length; k++){
                for(auto &item : aln_ret[k]){
                    int l = item.first;
                    if (tmp_aln_ret[k].find(l) == tmp_aln_ret[l].end()) continue; // N.B.
                    if (tmp_aln_ret[k][l].value >= 0.01){
                        new_aln_results[i][j][k][l].value = tmp_aln_ret[k][l].value;
                        new_aln_results[j][i][l][k].value = tmp_aln_ret[k][l].value;
                        // cout << i << " " << j << " " << k << " " << l << " " << tmp_aln_ret[k][l].value << endl;
                    }
                }
            }
            delete[] tmp_aln_ret;
        }
    }

    return new_aln_results;
}

MultiSequence *LinearAlignAlignments (MultiSequence *align1, MultiSequence *align2,
                                const vector<vector<unordered_map<int, multi_aln_ret>*>> &multi_aln_results,
                                const ProbabilisticModel &model, int hmmBeam){

    // Print some info about the alignment
    // SafeVector<float> *posterior = model.LinearBuildPosterior (align1, align2, multi_aln_results);
    // model.ComputeAlignment (align1->GetSequence(0)->GetLength(), align2->GetSequence(0)->GetLength(), *posterior);
    // delete posterior;

    // Choose the alignment routine depending on the "cosmetic" gap penalties used
    const unordered_map<int, multi_aln_ret> *sum_aln_ret = model.LinearMultiAlnResults(align1, align2, multi_aln_results);
    pair<SafeVector<char> *, float> alignment = model.LinearComputeAlignment(hmmBeam, align1->GetSequence(0)->GetLength(), align2->GetSequence(0)->GetLength(), sum_aln_ret);
    delete[] sum_aln_ret;

    // Build final alignment
    MultiSequence *result = new MultiSequence();
    for (int i = 0; i < align1->GetNumSequences(); i++)
        result->AddSequence (align1->GetSequence(i)->AddGaps(alignment.first, 'X'));
    for (int i = 0; i < align2->GetNumSequences(); i++)
        result->AddSequence (align2->GetSequence(i)->AddGaps(alignment.first, 'Y'));
    result->SortByLabel();

    // Free temporary alignment
    delete alignment.first;

    return result;
}

MultiSequence *LinearProcessTree (const TreeNode *tree, MultiSequence *sequences,
                            const vector<vector<unordered_map<int, multi_aln_ret>*>> &multi_aln_results,
                            const ProbabilisticModel &model, int hmmBeam){
    MultiSequence *result;

    // Check if this is a node of the alignment tree
    if (tree->GetSequenceLabel() == -1){
        MultiSequence *alignLeft = LinearProcessTree (tree->GetLeftChild(), sequences, multi_aln_results, model, hmmBeam);
        MultiSequence *alignRight = LinearProcessTree (tree->GetRightChild(), sequences, multi_aln_results, model, hmmBeam);

        assert (alignLeft);
        assert (alignRight);

        result = LinearAlignAlignments (alignLeft, alignRight, multi_aln_results, model, hmmBeam);
        assert (result);

        delete alignLeft;
        delete alignRight;
    }

    // Otherwise, this is a leaf of the alignment tree
    else {
        result = new MultiSequence(); assert (result);
        result->AddSequence (sequences->GetSequence(tree->GetSequenceLabel())->Clone());
    }

    return result;
}

void LinearDoIterativeRefinement (const vector<vector<unordered_map<int, multi_aln_ret>*>> &multi_aln_results,
                            const ProbabilisticModel &model, MultiSequence* &alignment, int i, int hmmBeam){
    set<int> groupOne, groupTwo;
    randomnumber rn;
    rn.seed(1234+i);
    // Create two separate groups
    for (int i = 0; i < alignment->GetNumSequences(); i++){
        int x = rn.roll_int(1,10);
        //cout << "rand: " << x << '\t';
        if (x % 2) {
          groupOne.insert (i);
        }
        else
          groupTwo.insert (i);
    }

    if (groupOne.empty() || groupTwo.empty()) return;

    // Project into the two groups
    MultiSequence *groupOneSeqs = alignment->Project (groupOne); assert (groupOneSeqs);
    MultiSequence *groupTwoSeqs = alignment->Project (groupTwo); assert (groupTwoSeqs);
    delete alignment;

    // Realign
    alignment = LinearAlignAlignments (groupOneSeqs, groupTwoSeqs, multi_aln_results, model, hmmBeam);

    delete groupOneSeqs;
    delete groupTwoSeqs;
}

MultiSequence *LinearComputeFinalAlignment (const TreeNode *tree, MultiSequence *sequences,
                                      const vector<vector<unordered_map<int, multi_aln_ret>*>> &multi_aln_results,
                                      const ProbabilisticModel &model, int hmmBeam){
    MultiSequence *alignment = LinearProcessTree (tree, sequences, multi_aln_results, model, hmmBeam);

    // Iterative refinement
    for (int i = 0; i < numIterativeRefinementReps; i++)
        LinearDoIterativeRefinement (multi_aln_results, model, alignment, i, hmmBeam);

    return alignment;
}