#include "ProbabilisticModel.h"
#include <list>
#include <cmath>
#include <cstdio>


using namespace std;

//! ProbabilisticModel Class
/*!
    The ProbabilisticModel Class stores the parameters of a probabilistic model.

*/


//! Store the largest of three values x1, x2, and x3 in *x.  
//! If xi is the largest value, then store bi in *b.
// void ChooseBestOfThree(float x1, float x2, float x3, char b1, char b2, char b3, float *x, char *b){
//     if (x1 >= x2){
//         if (x1 >= x3){
//             *x = x1;
//             *b = b1;
//             return;
//         }
//         *x = x3;
//         *b = b3;
//         return;
//     }
//     if (x2 >= x3){
//         *x = x2;
//         *b = b2;
//         return;
//     }
//     *x = x3;
//     *b = b3;
// }

unsigned long tmp_quickselect_partition(vector<pair<double, int>>& scores, unsigned long lower, unsigned long upper) {
    double pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}

double tmp_quickselect(vector<pair<double, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = tmp_quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return tmp_quickselect(scores, lower, split-1, k);
    else return tmp_quickselect(scores, split+1, upper, k - length);
}

double tmp_beam_prune(std::unordered_map<int, AlignState> &beamstep, int beamsize){
    vector<pair<double, int>> scores;
    for (auto &item : beamstep) {
        int ik = item.first;
        AlignState &cand = item.second;
        scores.push_back(make_pair(cand.alpha, ik));
        // scores.push_back(make_pair(cand.alpha, ik));
    }
    if (scores.size() <= beamsize) return VALUE_MIN;
    double threshold = tmp_quickselect(scores, 0, scores.size() - 1, scores.size() - beamsize);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }
    // cout << "threshold: " << threshold << endl;
    return threshold;
}

pair<SafeVector<char> *, float> ProbabilisticModel::LinearComputeAlignment(int hmmBeam, int seq1Length, int seq2Length, const unordered_map<int, multi_aln_ret>* &mul_aln_rets) const {
    unordered_map<int, multi_aln_ret>* tmp_mul_aln_rets = const_cast<unordered_map<int, multi_aln_ret>*>(mul_aln_rets);
    unordered_map<int, AlignState>* tmp_aln_rets = new unordered_map<int, AlignState>[seq1Length + seq2Length + 3];
    unsigned max_len = seq1Length > seq2Length ? seq1Length : seq2Length;
    // step 0
    tmp_aln_rets[0][0].alpha = 0;
    tmp_aln_rets[0][0].manner = 3; // ALIGN_ALN;
    tmp_aln_rets[0][0].step = 0;
    tmp_aln_rets[0][0].i = 0;
    tmp_aln_rets[0][0].k = 0;

    for(int s = 0; s < seq1Length + seq2Length + 1; ++s){
        int beamsize = hmmBeam;
        if (beamsize > 0 && tmp_aln_rets[s].size() > beamsize) tmp_beam_prune(tmp_aln_rets[s], beamsize);
        for (auto &item : tmp_aln_rets[s]) {
            AlignState &state = item.second;
            unsigned i = state.i;
            unsigned k = state.k;
            unsigned step = state.step;
            unsigned next_key;
            int manner = state.manner;
            int next_manner;
            unsigned next_i, next_k, next_step;
            for (int m = 3; m >= 1; m--){
                next_step = step;
                switch (m)
                {
                case 3: // ALIGN_ALN:
                    next_manner = m; // ALIGN_ALN;
                    next_i = i + 1;
                    next_k = k + 1;
                    next_step += 2;
                    next_key = next_i * max_len + next_k;
                    if ((next_i <= seq1Length) && (next_k <= seq2Length)) {
                        if (mul_aln_rets[next_i].find(next_k) == mul_aln_rets[next_i].end()) continue; // N.B.
                        tmp_aln_rets[next_step][next_key].alpha = tmp_mul_aln_rets[next_i][next_k].value; // N.B.
                        tmp_aln_rets[next_step][next_key].alpha += state.alpha; 
                        tmp_aln_rets[next_step][next_key].i = next_i;
                        tmp_aln_rets[next_step][next_key].k = next_k;
                        tmp_aln_rets[next_step][next_key].step = next_step;
                        tmp_aln_rets[next_step][next_key].manner = m;
                        // cout << i << " " << k << " " << next_i << " " << next_k << " " << tmp_aln_rets[next_step][next_key].alpha << endl;
                    }
                    break;

                case 1: // ALIGN_INS1:
                    next_manner = m; // ALIGN_INS1;
                    next_i = i + 1;
                    next_k = k;
                    next_step += 1;
                    next_key = next_i * max_len + next_k;
                    if ((next_i <= seq1Length) && (next_k <= seq2Length)){
                        if (state.alpha > tmp_aln_rets[next_step][next_key].alpha){
                            tmp_aln_rets[next_step][next_key].alpha = state.alpha; // N.B.
                            tmp_aln_rets[next_step][next_key].i = next_i;
                            tmp_aln_rets[next_step][next_key].k = next_k;
                            tmp_aln_rets[next_step][next_key].step = next_step;
                            tmp_aln_rets[next_step][next_key].manner = m;
                            // cout << next_i << " " << next_k << " " << tmp_aln_rets[next_step][next_key].alpha << endl;
                        }
                    }
                    break;

                case 2: // ALIGN_INS2:
                    next_manner = m; // ALIGN_INS2;
                    next_i = i;
                    next_k = k + 1;
                    next_step += 1;
                    next_key = next_i * max_len + next_k;
                    if ((next_i <= seq1Length) && (next_k <= seq2Length)) {
                        if (state.alpha > tmp_aln_rets[next_step][next_key].alpha){
                            tmp_aln_rets[next_step][next_key].alpha = state.alpha; // N.B.
                            tmp_aln_rets[next_step][next_key].i = next_i;
                            tmp_aln_rets[next_step][next_key].k = next_k;
                            tmp_aln_rets[next_step][next_key].step = next_step;
                            tmp_aln_rets[next_step][next_key].manner = m;
                            // cout << next_i << " " << next_k << " " << tmp_aln_rets[next_step][next_key].alpha << endl;
                        }
                    }
                    break;
                default:
                    break;
                }
            }
        }
    }

    // compute traceback
    SafeVector<char> *alignment = new SafeVector<char>; 
    assert (alignment);

    int seq1_pos = seq1Length;
    int seq2_pos = seq2Length;
    int step = seq1_pos + seq2_pos;
    int key = seq1_pos * max_len + seq2_pos;
    int pre_manner = tmp_aln_rets[seq1Length + seq2Length][key].manner;
    double total = tmp_aln_rets[seq1Length + seq2Length][key].alpha;
    // cout << "tmp_aln_rets[seq1Length + seq2Length][key].manner" << pre_manner << endl;
    while (true) {
        if ((seq1_pos == 0) && (seq2_pos == 0)) break;
        // cout << step << " " << key << " " << seq1_pos << " " << seq2_pos << " " << pre_manner  << endl;
        switch (pre_manner)
        {
        case 3:
            alignment->push_back ('B');
            // cout << "B";
            seq1_pos--;
            seq2_pos--;
            step -= 2;
            key = seq1_pos * max_len + seq2_pos;
            pre_manner = tmp_aln_rets[step][key].manner;
            break;
        case 2: // INS2
            alignment->push_back ('Y');
            // cout << "Y";
            seq2_pos--;
            step--;
            key = seq1_pos * max_len + seq2_pos;
            pre_manner = tmp_aln_rets[step][key].manner;
            break;
        case 1: // INS1
            alignment->push_back ('X');
            // cout << "X";
            seq1_pos--;
            step--;
            key = seq1_pos * max_len + seq2_pos;
            pre_manner = tmp_aln_rets[step][key].manner;
            break;
        case 0:
            cout << "ComputeAlignment traceback error !!" << endl;
            cout << seq1_pos << " " << seq2_pos << endl;
            exit(0);
        default:
            break;
        }
    }
    // cout << endl;
    // alignment->push_back('X');
    reverse(alignment->begin(), alignment->end());
    // for(auto i = alignment->begin(); i!=alignment->end(); i++){
    //     cout << *i;
    // }
    // cout << endl;

    // cout << "total: " << total - mul_aln_rets[0][0].value << endl;
    delete[] tmp_aln_rets;
    return make_pair(alignment, total);
}

unordered_map<int, multi_aln_ret> * ProbabilisticModel::LinearMultiAlnResults(MultiSequence *align1, MultiSequence *align2, const vector<vector<unordered_map<int, multi_aln_ret>*>> &mul_aln_rets, float cutoff) const {
    const int seq1Length = align1->GetSequence(0)->GetLength();
    const int seq2Length = align2->GetSequence(0)->GetLength();

    unordered_map<int, multi_aln_ret>* sum_aln_ret = new unordered_map<int, multi_aln_ret>[seq1Length + 1];
    // cout << seq1Length << " " << seq2Length << endl;
    for (int i = 0; i < align1->GetNumSequences(); i++){
        int first = align1->GetSequence(i)->GetLabel();
        SafeVector<int> *mapping1 = align1->GetSequence(i)->GetMapping();
        // Loops through align2
        for (int j = 0; j < align2->GetNumSequences(); j++){
            int second = align2->GetSequence(j)->GetLabel();
            SafeVector<int> *mapping2 = align2->GetSequence(j)->GetMapping();
            // cout << "seqs: " << i << " " << j << " " << align1->GetSequence(i)->GetLength() << " " << align2->GetSequence(j)->GetLength() << endl;

            if(first < second){
                unordered_map<int, multi_aln_ret>* aln_ret = mul_aln_rets[first][second];
                int seq1len = mapping1->size() - 1;
                for (int ii = 0; ii <= seq1len; ii++){
                    int ibase = (*mapping1)[ii];
                    for (auto &item : aln_ret[ii]) {
                        // cout << "ii: " << ii << " " << item.first << endl;
                        int jbase = (*mapping2)[item.first];
                        if (item.second.value < 0.01) continue;
                        sum_aln_ret[ibase][jbase].value += item.second.value;
                        // cout << i <<  " "  << j <<  " "  << first << " "  << second << " " << ii << " "  << item.first << " " << ibase  << " " << jbase  << " " <<  item.second.value << " " << sum_aln_ret[ibase][jbase].value << endl;
                    }
                }
            } else {
                unordered_map<int, multi_aln_ret>* aln_ret = mul_aln_rets[second][first];
                int seq2len = mapping2->size() - 1;
                for (int jj = 0; jj <= seq2len; jj++){
                    int jbase = (*mapping2)[jj];
                    for (auto &item : aln_ret[jj]) {
                        int ibase = (*mapping1)[item.first];
                        if (item.second.value < 0.01) continue;
                        sum_aln_ret[ibase][jbase].value += item.second.value;
                        // cout << i <<  " "  << j <<  " "  <<  first << " "  << second << " " << ibase  << " " << jbase <<  item.second.value << " " << sum_aln_ret[ibase][jbase].value << endl;
                    }
                }
            }
            delete mapping2;
        }
        delete mapping1;
    }

    return sum_aln_ret;
}

SafeVector<float> * ProbabilisticModel::LinearBuildPosterior(MultiSequence *align1, MultiSequence *align2,
                     const vector<vector<unordered_map<int, multi_aln_ret>*>> &mul_aln_rets, float cutoff) const {
    const int seq1Length = align1->GetSequence(0)->GetLength();
    const int seq2Length = align2->GetSequence(0)->GetLength();
    // cout << "seq1Length: " << seq1Length << ", seq2Length: " << seq2Length << endl;
    // cout << "align1 size: " << align1->GetNumSequences() << ", align2 size: " << align2->GetNumSequences()  << endl;

    SafeVector<float> *posteriorPtr = new SafeVector<float>((seq1Length+1) * (seq2Length+1), 0); assert (posteriorPtr);
    SafeVector<float> &posterior = *posteriorPtr;
    SafeVector<float>::iterator postPtr = posterior.begin();

    for (int i = 0; i < align1->GetNumSequences(); i++){
        int first = align1->GetSequence(i)->GetLabel();
        SafeVector<int> *mapping1 = align1->GetSequence(i)->GetMapping();

        // Loops through align2
        for (int j = 0; j < align2->GetNumSequences(); j++){
            int second = align2->GetSequence(j)->GetLabel();
            SafeVector<int> *mapping2 = align2->GetSequence(j)->GetMapping();
            // cout << "first: " << first << ", second: " << second << endl;

            if(first < second){
                // get the associated alignment posterior probs.
                unordered_map<int, multi_aln_ret>* aln_ret = mul_aln_rets[first][second];
                
                int seq1len = mapping1->size() - 1;
                // cout << "seq1 length: " << seq1len << endl;
                for (int ii = 1; ii <= seq1len; ii++){
                    int base = (*mapping1)[ii] * (seq2Length+1);
                    // add in all relevant values
                    for (auto &item : aln_ret[ii]) {
                        int jj = item.first;
                        // cout << ii << " " << base << " " << jj << endl;
                        posterior[base + (*mapping2)[jj]] += item.second.value;
                        // cout << ii << " " << base << " " << jj << " " << (*mapping2)[jj] <<  " " << base + (*mapping2)[jj] << " " << item.second.value << endl;
                    }
                }
            } else {
                // get the associated alignment posterior probs.
                unordered_map<int, multi_aln_ret>* aln_ret = mul_aln_rets[second][first];
                int seq2len = mapping2->size() - 1;
                // cout << "seq2 length: " << seq2len << endl;
                for (int jj = 1; jj <= seq2len; jj++){
                    int base = (*mapping2)[jj];
                    // add in all relevant values
                    for (auto &item : aln_ret[jj]) {
                        int ii = item.first;
                        // cout << jj << " " << ii << endl;
                        posterior[base + (*mapping1)[ii] * (seq2Length + 1)] += item.second.value;
                        // cout << jj << " " << base << " " << ii << " " << (*mapping1)[ii] <<  " " << (*mapping1)[ii] * (seq2Length + 1) << " " <<  base + (*mapping1)[ii] * (seq2Length + 1) << " " << item.second.value << endl;
                    }
                }
            }
            delete mapping2;
        }
        delete mapping1;
    }
    return posteriorPtr;
}
