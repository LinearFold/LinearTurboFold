/*
 *LinearAlignment.cpp*
 The main code for LinearAlignment: Linear-Time Approximation of 
                                    pairwise RNA sequences alignment
                                    and alignment co-incidence probabilities

 author: Sizhen Li
 created by: 06/2020
*/

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string.h>
#include <map>
#include <stdio.h> 
#include <set> 

#include "LinearAlign.h"
#include "../../LinearTurboFold.h"

using namespace std;

int partition(vector<int>& A, int p,int q)
{
    int x= A[p];
    int i=p;
    int j;

    for(j=p+1; j<q; j++)
    {
        if(A[j]<=x)
        {
            i=i+1;
            swap(A[i],A[j]);
        }
    }
    swap(A[i],A[p]);
    return i;
}

void quickSort(vector<int>& A, int p,int q)
{
    int r;
    if(p<q)
    {
        r=partition(A, p,q);
        quickSort(A,p,r);  
        quickSort(A,r+1,q);
    }
}

unsigned long BeamAlign::quickselect_partition(vector<pair<double, int>>& scores, unsigned long lower, unsigned long upper) {
    double pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}

// in-place quick-select
double BeamAlign::quickselect(vector<pair<double, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}

double BeamAlign::beam_prune(std::unordered_map<int, AlignState> &beamstep){
    scores.clear();
    for (auto &item : beamstep) {
        int ik = item.first;
        AlignState &cand = item.second;
        scores.push_back(make_pair(cand.alpha, ik));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    double threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}

void update_if_better(AlignState &state, double newscore, int pre_manner, int manner, unsigned step, unsigned i, unsigned k) {
    if (state.alpha < newscore)
        // ++ nos_set_update;
        state.set(newscore, pre_manner, manner, step, i, k);
};

void update(AlignState &state, double newscore, int pre_manner, int manner, unsigned step, unsigned i, unsigned k) {
    state.set(newscore, pre_manner, manner, step, i, k);
};

double BeamAlign::get_trans_emit_prob(int prev_state, int current_state, int i, int k, double** &trans_probs, double** &emit_probs){
    // get_trans_prob
    if (prev_state < 0 || current_state < 0) return xlog(0.0); // N.B.
    double trans_prob = trans_probs[prev_state][current_state];

    // get_emit_prob
    int i_sym;
	int k_sym;

	// Fix symbols to gaps in case of of insertions.
	if(current_state == 0 || k == 0)
	{
		// Gap is coded into value 4 in the emission table.
		k_sym = 4; 
	}
	else
	{
		k_sym = nucs2[k];
	}

	if(current_state == 1 || i == 0)
	{
		// Gap is coded into value 4 in the emission table.
		i_sym = 4;
	}
	else
	{
		i_sym = nucs1[i];
	}

	// Compute the symbol index into emission table using the coded nucleotide values:
	// A->0, C->1, G->2, U->3, T->3, .->4
	// This defines a counting system in base of 5. (25 values.)
	// There are also emission of start and end symbols. These correspond to 25th and 26th indices in the emission probability table.
	int sym_index = i_sym * 5 + k_sym;

	// Check for exceptional cases of start and end symbols.
	// The indices correspond to the start symbol?
	if(i == 0 && k == 0)
	{
		sym_index = 25;
	}

	// The indices correspond to the end symbol?
	if(i == seq1_len && k == seq2_len)
	{
		sym_index = 26;
	}
    double emit_prob = emit_probs[sym_index][current_state];

	return(xlog_mul(emit_prob, trans_prob));
}

double BeamAlign::get_match_prior(int i, int k, bool prior)
{
    if (prior){
        if(i == 0 || k == 0 || i == seq1_len || k == seq2_len)
            return(0.0f); 
        
        return(xlog(out->get_match_score(i_seq1, i_seq2, i, k)));
    }
    
    return 0.0f;
}

void BeamAlign::prepare(string &seq1, string &seq2) {
    seq1_len = static_cast<unsigned>(seq1.length() + 1);
    seq2_len = static_cast<unsigned>(seq2.length() + 1);
    max_len = seq1_len >= seq2_len ? seq1_len : seq2_len;

    nucs1 = new int[seq1_len];
    nucs2 = new int[seq2_len];

    bestALN = new unordered_map<int, AlignState>[seq1_len + seq2_len + 1];
    bestINS1 = new unordered_map<int, AlignState>[seq1_len + seq2_len + 1];
    bestINS2 = new unordered_map<int, AlignState>[seq1_len + seq2_len + 1];

    scores.reserve(seq2_len);

    replace(seq1.begin(), seq1.end(), 'T', 'U');
    replace(seq2.begin(), seq2.end(), 'T', 'U');
    for (int i = 0; i < seq1_len; ++i){
        if (i == 0) nucs1[i] = 0;
        else nucs1[i] = GET_ACGU_NUM(seq1[i-1]);
        if (nucs1[i] > 3) {
            int r = rand() % 4;
            nucs1[i] = r;
        }
    }
    
    for (int i = 0; i < seq2_len; ++i){
        if (i == 0) nucs2[i] = 0;
        else nucs2[i] = GET_ACGU_NUM(seq2[i-1]);
        if (nucs2[i] > 3) {
            int r = rand() % 4;
            nucs2[i] = r;
        }
    }
}

void BeamAlign::traceback(vector<char> &aln1, vector<char> &aln2){
    unsigned i_seq1 = seq1_len;
    unsigned i_seq2 = seq2_len;
    unsigned step = i_seq1 + i_seq2;
    unsigned key = i_seq1 * max_len + i_seq2;
    int cur_manner = bestALN[step][key].manner;
    double best_score =  bestALN[step][key].alpha;

    while (true) {
        if ((i_seq1 == 0) && (i_seq2 == 0)) break;
        switch (cur_manner) {
            case 3: // ALIGN_ALN:
                if (i_seq1 < seq1_len) aln1.push_back(GET_NUC(nucs1[i_seq1]));
                if (i_seq2 < seq2_len) aln2.push_back(GET_NUC(nucs2[i_seq2]));
                cur_manner = bestALN[step][key].pre;
                step = step - 2;
                i_seq1 --;
                i_seq2 --;
                key = i_seq1 * max_len + i_seq2;
                break;
            case 1: // ALIGN_INS1:
                aln1.push_back(GET_NUC(nucs1[i_seq1]));
                aln2.push_back('-');
                cur_manner = bestINS1[step][key].pre;
                step = step - 1;
                i_seq1 --;
                key = i_seq1 * max_len + i_seq2;
                break;
            case 2: //ALIGN_INS2:
                aln1.push_back('-');
                aln2.push_back(GET_NUC(nucs2[i_seq2]));
                cur_manner = bestINS2[step][key].pre;
                step = step - 1;
                i_seq2 --;
                key = i_seq1 * max_len + i_seq2;
                break;
            default:  // MANNER_NONE or other cases
                printf("wrong manner at %d, %d: manner %d\n", i_seq1, i_seq2, cur_manner); fflush(stdout);
                assert(false);
        }
    }

    // std::reverse(aln1.begin(), aln1.end());
    // std::reverse(aln2.begin(), aln2.end());
    // for(auto e : aln1) cout << e;
    // cout << endl;
    // for(auto e : aln2) cout << e;
    // cout << endl;

    // return get_aln_similarity(aln1, aln2, '-');
}

void BeamAlign::ml_alignment(string &seq1, string &seq2, vector<char> &aln1, vector<char> &aln2, double** &transprobs, double** &emitprobs, bool prior){
    prepare(seq1, seq2);
        
    double trans_emit_prob;
    for(int s = 0; s < seq1_len + seq2_len; ++s){
        int nuc1 = nucs1[s];
        unordered_map<int, AlignState>& beamALN = bestALN[s];
        unordered_map<int, AlignState>& beamINS1 = bestINS1[s];
        unordered_map<int, AlignState>& beamINS2 = bestINS2[s];

        // initial state
        if (s == 0){
            beamALN[0].alpha = xlog(1.0);
            beamALN[0].manner = 3; // ALIGN_ALN;
            beamALN[0].step = 0;
            beamALN[0].i = 0;
            beamALN[0].k = 0;
        }

        vector<unordered_map<int, AlignState>*> beams{&beamINS1, &beamINS2, &beamALN};
        for (int i=0; i < beams.size(); i++){
            unordered_map<int, AlignState> &beamstep = *beams[i];

            if (beam > 0 && beamstep.size() > beam) beam_prune(beamstep);
            for (auto &item : beamstep) {
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
                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))) {
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                            trans_emit_prob = xlog_mul(get_match_prior(next_i, next_k, prior), trans_emit_prob); // match score
                            update_if_better(bestALN[next_step][next_key], xlog_mul(state.alpha, trans_emit_prob), manner, next_manner, next_step, next_i, next_k);
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = m; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step += 1;
                        next_key = next_i * max_len + next_k;
                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))) {
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                            update_if_better(bestINS1[next_step][next_key], xlog_mul(state.alpha, trans_emit_prob), manner, next_manner, next_step, next_i, next_k);
                        }
                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = m; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step += 1;
                        next_key = next_i * max_len + next_k;
                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))) {
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                            update_if_better(bestINS2[next_step][next_key], xlog_mul(state.alpha, trans_emit_prob), manner, next_manner, next_step, next_i, next_k);
                        }
                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }
    }

    double forward_score = bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].alpha;
    traceback(aln1, aln2);
}


double BeamAlign::forward(string seq1, string seq2, double** &trans_probs, double** &emit_probs, bool prior){
    // re-new 
    delete[] bestINS1;
    delete[] bestINS2;
    delete[] bestALN;

    bestALN = new unordered_map<int, AlignState>[seq1_len + seq2_len + 1];
    bestINS1 = new unordered_map<int, AlignState>[seq1_len + seq2_len + 1];
    bestINS2 = new unordered_map<int, AlignState>[seq1_len + seq2_len + 1];

    scores.reserve(seq2_len);

    // initial state
    bestALN[0][0].alpha = xlog(1.0);
    bestALN[0][0].manner = 3; // ALIGN_ALN;
    bestALN[0][0].step = 0;
    bestALN[0][0].i = 0;
    bestALN[0][0].k = 0;

    double trans_emit_prob;
    for(int s = 0; s < seq1_len + seq2_len; ++s){
        unordered_map<int, AlignState> &beamALN = bestALN[s];
        unordered_map<int, AlignState> &beamINS1 = bestINS1[s];
        unordered_map<int, AlignState> &beamINS2 = bestINS2[s];

        vector<unordered_map<int, AlignState>*> beams{&beamINS1, &beamINS2, &beamALN};

        for (int i=0; i < beams.size(); i++){
            unordered_map<int, AlignState> &beamstep = *beams[i];
            if (beam > 0 && beamstep.size() > beam) beam_prune(beamstep);
            for (auto &item : beamstep) {
                AlignState &state = item.second;
                unsigned i = state.i;
                unsigned k = state.k;
                unsigned step = state.step;
                int manner = state.manner;

                int next_manner;
                unsigned next_i, next_k, next_step, next_key;
                for (int m = 1; m <= 3; m++){
                    next_step = step;
                    switch (m)
                    {
                    case 3: // ALIGN_ALN:
                        next_manner = m; // ALIGN_ALN;
                        next_i = i + 1;
                        next_k = k + 1;
                        next_step += 2;
                        
                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))) {
                            next_key = next_i * max_len + next_k;
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, trans_probs, emit_probs);
                            trans_emit_prob = xlog_mul(get_match_prior(next_i, next_k, prior), trans_emit_prob);
                            double newscore = xlog_sum(bestALN[next_step][next_key].alpha, xlog_mul(state.alpha, trans_emit_prob));
                            update(bestALN[next_step][next_key], newscore, manner, next_manner, next_step, next_i, next_k);
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = m; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step += 1;

                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))){
                            next_key = next_i * max_len + next_k;
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, trans_probs, emit_probs);
                            double newscore = xlog_sum(bestINS1[next_step][next_key].alpha, xlog_mul(state.alpha, trans_emit_prob));
                            update(bestINS1[next_step][next_key], newscore, manner, next_manner, next_step, next_i, next_k);
                        }
                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = m; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step += 1;

                        if (((next_i < seq1_len) && (next_k < seq2_len)) || ((next_i == seq1_len) && (next_k == seq2_len))) {
                            next_key = next_i * max_len + next_k;
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, trans_probs, emit_probs);
                            double newscore = xlog_sum(bestINS2[next_step][next_key].alpha, xlog_mul(state.alpha, trans_emit_prob));
                            update(bestINS2[next_step][next_key], newscore, manner, next_manner, next_step, next_i, next_k);
                        }
                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }
    }

    double forward_score = bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].alpha;
    // printf("forward score: %.5f\n", forward_score);
    return forward_score;
}

double BeamAlign::backward(double** &transprobs, double** &emitprobs, bool prior){
    bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].beta = xlog(1.0);
    bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].manner = 3; //ALIGN_ALN;
    bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].step = seq1_len + seq2_len;
    bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].i = seq1_len;
    bestALN[seq1_len + seq2_len][seq1_len * max_len + seq2_len].k = seq2_len;

    for(int s = seq1_len + seq2_len - 2; s >= 0; --s) {
        double trans_emit_prob;
        unordered_map<int, AlignState> &beamALN = bestALN[s];
        unordered_map<int, AlignState> &beamINS1 = bestINS1[s];
        unordered_map<int, AlignState> &beamINS2 = bestINS2[s];
        vector<unordered_map<int, AlignState>*> beams{&beamALN, &beamINS1, &beamINS2}; // N.B.

        for (int i=0; i < beams.size(); i++){
            unordered_map<int, AlignState> &beamstep = *beams[i]; // N.B.
            for (auto &item : beamstep) {
                AlignState &state = item.second;
                unsigned i = state.i;
                unsigned k = state.k;
                unsigned step = state.step;
                int manner = state.manner;

                int next_manner;
                unsigned next_i, next_k, next_step, next_key;
                for (int m = 1; m <= 3; m++){
                    next_step = step;
                    switch (m)
                    {
                    case 3: // ALIGN_ALN:
                        next_manner = m; // ALIGN_ALN;
                        next_i = i + 1;
                        next_k = k + 1;
                        next_step += 2;

                        if (((next_i == seq1_len) && (next_k == seq2_len)) || ((next_i < seq1_len) && (next_k < seq2_len))) {
                            next_key = next_i * max_len + next_k;
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                            trans_emit_prob = xlog_mul(get_match_prior(next_i, next_k, prior), trans_emit_prob);
                            state.beta = xlog_sum(state.beta, xlog_mul(bestALN[next_step][next_key].beta, trans_emit_prob));
                        }
                        break;

                    case 1: // ALIGN_INS1:
                        next_manner = m; // ALIGN_INS1;
                        next_i = i + 1;
                        next_k = k;
                        next_step += 1;

                        if ((next_i < seq1_len) && (next_k < seq2_len)) {
                            next_key = next_i * max_len + next_k;
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                            state.beta = xlog_sum(state.beta, xlog_mul(bestINS1[next_step][next_key].beta, trans_emit_prob));
                        }
                        break;

                    case 2: // ALIGN_INS2:
                        next_manner = m; // ALIGN_INS2;
                        next_i = i;
                        next_k = k + 1;
                        next_step += 1;

                        if ((next_i < seq1_len) && (next_k < seq2_len)) {
                            next_key = next_i * max_len + next_k;
                            trans_emit_prob = get_trans_emit_prob(manner - 1, next_manner - 1, next_i, next_k, transprobs, emitprobs);
                            state.beta = xlog_sum(state.beta, xlog_mul(bestINS2[next_step][next_key].beta, trans_emit_prob));
                        }
                        break;
                    
                    default:
                        break;
                    }
                }
            }
        }
    }

    double back_score = bestALN[0][0].beta;
    // printf("backward score: %.5f\n", back_score);
    return back_score;
}

std::unordered_map<int, aln_ret>*  BeamAlign::cal_align_prob(double forward_score, double threshold, std::unordered_map<int, aln_ret>* &aln_results){
    double aln_prob, ins1_prob, ins2_prob;
   
    for (int s = 0; s < seq1_len + seq2_len; s++){
        for(auto &item : bestALN[s]){
            AlignState &state = item.second;
            int i = state.i;
            int k = state.k;
            aln_prob = xlog_div(xlog_mul(state.alpha, state.beta), forward_score);
            if (aln_prob > float(-9.91152)) {
                aln_results[i][k].prob = xlog_sum(aln_results[i][k].prob, aln_prob);
                aln_results[i][k].aln_prob = aln_prob;
            }
        }

        for(auto &item : bestINS1[s]){
            AlignState &state = item.second;
            int i = state.i;
            int k = state.k;
            ins1_prob = xlog_div(xlog_mul(state.alpha, state.beta), forward_score);
            if (ins1_prob > float(-9.91152)) aln_results[i][k].prob = xlog_sum(aln_results[i][k].prob, ins1_prob);
        }

        for(auto &item : bestINS2[s]){
            AlignState &state = item.second;
            int i = state.i;
            int k = state.k;
            ins2_prob = xlog_div(xlog_mul(state.alpha, state.beta), forward_score);
            if (ins2_prob > float(-9.91152)) aln_results[i][k].prob = xlog_sum(aln_results[i][k].prob, ins2_prob);
        }
    }

    vector<int> cands;
    for (int i = 1; i < seq1_len; i++) {
        cands.clear();
        for(auto &item : aln_results[i]){
            int k = item.first;
            double prob = item.second.prob;
            if (prob < threshold) {
                cands.push_back(k);
            } else {
                item.second.prob = Fast_Exp(prob);
                item.second.aln_prob = Fast_Exp(item.second.aln_prob);
            }
        }
     
        for (auto &k : cands) {
            aln_results[i].erase(k);
        }
    }
    return aln_results;
}

BeamAlign::BeamAlign(int beam_size)
    : beam(beam_size){
}

BeamAlign::~BeamAlign(){
    delete[] bestINS1;
    delete[] bestINS2;
    delete[] bestALN;

    delete[] nucs1;
    delete[] nucs2;
}