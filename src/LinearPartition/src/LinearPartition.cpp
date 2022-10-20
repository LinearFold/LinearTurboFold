/*
 *LinearPartition.cpp*
 The main code for LinearPartition: Linear-Time Approximation of 
                                    RNA Folding Partition Function 
                                    and Base Pairing Probabilities

 author: He Zhang
 created by: 03/2019
*/

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <stdio.h> 
#include <set> 

#include "LinearPartition.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"
#include "bpp.h"
#include "../../SeqFold.h"
#include "../../utils/pfunction_math.h"
#include "../../utils/math/xlog_math.h"
#include "../../LinearTurboFold.h"
#define SPECIAL_HP


#define pf
// #define lpv // vinnea

using namespace std;

unsigned long quickselect_partition(vector<pair<float, int>>& scores, unsigned long lower, unsigned long upper) {
    float pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}

// in-place quick-select
float quickselect(vector<pair<float, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}


float BeamCKYParser::beam_prune(std::unordered_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        float newalpha = (k >= 0 ? bestC[k].alpha : 0.0) + cand.alpha;
        scores.push_back(make_pair(newalpha, i));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    float threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}

float BeamCKYParser::beam_prune(int i_seq, LinearTurboFold* outer, std::unordered_map<int, State> &beamstep, int j, unordered_map<int, ExtValue>* &ext_info) {
    if (beam < 0 || beamstep.size() <= beam) { 
        auto it = beamstep.begin();
        while (it != beamstep.end())
        {
            int i = it->first;
            double ex_info = outer->get_folding_extrinsic_information(i_seq, j+1, i+1);
            if (ex_info <= LOG_OF_ZERO) {
                it = beamstep.erase(it);
            } else {
                ext_info[j+1][i+1].score = 0.3 * ex_info;
                beamstep[i].alpha = PROD(beamstep[i].alpha, ext_info[j+1][i+1].score);
                it ++;
            } 
        }
    } else {
        scores.clear();
        for (auto &item : beamstep) {
            int i = item.first;
            State &cand = item.second;
            int k = i - 1;

            double ex_info = outer->get_folding_extrinsic_information(i_seq, j+1, i+1);
            if (ex_info <= LOG_OF_ZERO) {
                ext_info[j+1][i+1].score = LOG_OF_ZERO;
                cand.alpha = LOG_OF_ZERO;
            } else {
                ext_info[j+1][i+1].score = 0.3 * ex_info;
                cand.alpha = PROD(cand.alpha, ext_info[j+1][i+1].score);
            }
            float newalpha = PROD((k >= 0 ? bestC[k].alpha : 0.0), (ex_info>LOG_OF_ZERO) ? cand.alpha : LOG_OF_ZERO);
            scores.push_back(make_pair(newalpha, i));
        }
        
        float threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
        for (auto &p : scores) {
            if ((p.first < threshold) || (p.first <= -700000.0f)) { // lisiz: log_of_zero (-709783)
                beamstep.erase(p.second);
            }
        } 
    }

    return VALUE_MIN;
}

void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;

    if (allocate) {
        delete[] nucs;
        delete[] bestC;
        delete[] bestH;
        // delete[] bestP;
        delete[] bestM;
        delete[] bestM2;
        delete[] bestMulti;
    }
    allocate = true;

    nucs = new int[seq_length];
    bestC = new State[seq_length];
    bestH = new unordered_map<int, State>[seq_length];
    // bestP = new unordered_map<int, State>[seq_length];
    bestM = new unordered_map<int, State>[seq_length];
    bestM2 = new unordered_map<int, State>[seq_length];
    bestMulti = new unordered_map<int, State>[seq_length];
    
    scores.reserve(seq_length);
}


double BeamCKYParser::parse(int i_iter, int i_seq, LinearTurboFold* outer, string& seq, unordered_map<int, State>* &pfscore, unordered_map<int, ExtValue>* &ext_info, string pfSaveFile, string bppSaveFile) {
    // convert T to U
    replace(seq.begin(), seq.end(), 'T', 'U');

    prepare(static_cast<unsigned>(seq.length()));

    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    vector<int> next_pair[NOTON];
    {
        for (int nuci = 0; nuci < NOTON; ++nuci) {
            // next_pair
            next_pair[nuci].resize(seq_length, -1);
            int next = -1;
            for (int j = seq_length-1; j >=0; --j) {
                next_pair[nuci][j] = next;
                if (_allowed_pairs[nuci][nucs[j]]) next = j;
            }
        }
    }


#ifdef SPECIAL_HP
    v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#endif
    if(seq_length > 0) bestC[0].alpha = 0.0;
    if(seq_length > 1) bestC[1].alpha = 0.0;

    value_type newscore;
    for(int j = 0; j < seq_length; ++j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = pfscore[j]; // bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of H
        {
            if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);

            {
                // for nucj put H(j, j_next) into H[j_next]
                int jnext = next_pair[nucj][j];
                if (no_sharp_turn) while (jnext - j < 4 && jnext != -1) jnext = next_pair[nucj][jnext];
                if (jnext != -1) {
                    int nucjnext = nucs[jnext];
                    int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-j-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[j];
                        else if (jnext-j-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[j];
                        else if (jnext-j-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[j];
#endif
                        newscore = - v_score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
                        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore/kT);
                }
            }

            {
                // for every state h in H[j]
                //   1. extend h(i, j) to h(i, jnext)
                //   2. generate p(i, j)
                for (auto &item : beamstepH) {
                    int i = item.first;
                    State &state = item.second;
                    int nuci = nucs[i];
                    int jnext = next_pair[nuci][j];

                    if (jnext != -1) {
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucjnext = nucs[jnext];
                        int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                        // 1. extend h(i, j) to h(i, jnext)
                        // value_type newscore;

                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-i-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[i];
                        else if (jnext-i-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[i];
                        else if (jnext-i-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[i];
#endif
                        newscore = - v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
                        Fast_LogPlusEquals(bestH[jnext][i].alpha, (newscore/kT));
                    }

                    // 2. generate p(i, j)
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha);  
                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);

            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
                        Fast_LogPlusEquals(bestMulti[jnext][i].alpha, (state.alpha));
                    }
                }

                // 2. generate P (i, j)
                {
                    value_type score_multi = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(beamstepP[i].alpha, (state.alpha + score_multi/kT));
                }
            }
        }

        // beam of P
        {   
            if (i_iter == 0){
                if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);
            } else {
                beam_prune(i_seq, outer, beamstepP, j, ext_info);
            }

            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;
                
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 1. generate new helix / single_branch
                // new state is of shape p..i..j..q
                if (i >0 && j<seq_length-1) {
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
                                int score_single = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                            nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(pfscore[q][p].alpha, (state.alpha + score_single/kT));
                                // Fast_LogPlusEquals(bestP[q][p].alpha, (state.alpha + score_single/kT));
                            } else {
                                // single branch
                                int score_single = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(pfscore[q][p].alpha, (state.alpha + score_single/kT));
                                // Fast_LogPlusEquals(bestP[q][p].alpha, (state.alpha + score_single/kT));
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }

                // 2. M = P
                if(i > 0 && j < seq_length-1){
                        int score_M1 = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepM[i].alpha, (state.alpha + score_M1/kT));
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
                    int M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    float m1_alpha = state.alpha + M1_score/kT;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(beamstepM2[newi].alpha, m_state.alpha + m1_alpha);
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                      State& prefix_C = bestC[k];
                        int nuck = nuci_1;
                        int nuck1 = nuci;
                        int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                             nucj, nucj1, seq_length);
                        state.score_external = score_external_paired/kT + prefix_C.alpha;                   
                        Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + score_external_paired/kT);      
                    } else {
                        int score_external_paired = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length);
                        state.score_external = score_external_paired/kT;
                        Fast_LogPlusEquals(beamstepC.alpha, state.alpha + score_external_paired/kT); 
                    }
                }
            }
        }

        // beam of M2
        {
            if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];
                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
                        Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha);
                        }
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha);  
            }
        }

        // beam of M
        {
            // float threshold = VALUE_MIN;
            // if (beam > 0 && beamstepM.size() > beam) threshold = beam_prune(beamstepM);
            if (beam > 0 && beamstepM.size() > beam) beam_prune(beamstepM);

            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha); 
                }
            }
        }

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
                Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha); 
            }
        }

    }  // end of for-loo j
    
    State& viterbi = bestC[seq_length-1];

    if(is_verbose) printf("Iter %d Seq %d Free Energy of Ensemble: %.10f kcal/mol\n", i_iter, i_seq, -kT * viterbi.alpha / 100.0);

    fflush(stdout);

    // outside phrase
    outside(outer, i_seq, i_iter, next_pair, pfscore, ext_info);

    if (!bppSaveFile.empty()){ 
        bpp_file = bppSaveFile;
        cal_PairProb(viterbi, pfscore, ext_info);
    }

    if (!pfSaveFile.empty()){ 
        forest_file = pfSaveFile;
        dump_forest(seq, pfscore, true);
    }

    return (double)viterbi.alpha;
}

// dump forest
void BeamCKYParser::print_states(FILE *fptr, unordered_map<int, State>& states, int j, string label, bool inside_only, double threshold) {    
    for (auto & item : states) {
        int i = item.first;
        State & state = item.second;
        if (inside_only)
            fprintf(fptr, "%s %d %d %.5lf\n", label.c_str(), i+1, j+1, state.alpha);
        else
            if (state.alpha + state.beta > threshold) // lhuang : alpha + beta - totalZ < ...
                fprintf(fptr, "%s %d %d %.5lf %.5lf\n", label.c_str(), i+1, j+1, state.alpha, state.beta);
    }
}

// dump forest
void BeamCKYParser::dump_forest(string seq, unordered_map<int, State>* &pfscore, bool inside_only) {  
//   printf("Dumping (%s) Forest to %s...\n", (inside_only ? "Inside-Only" : "Inside-Outside"), forest_file.c_str());
  FILE *fptr = fopen(forest_file.c_str(), "w");  // lhuang: should be fout >>
  fprintf(fptr, "%s\n", seq.c_str());
  int n = seq.length(), j;
  for (j = 0; j < n; j++) {
    if (inside_only)
        fprintf(fptr, "E %d %.5lf\n", j+1, bestC[j].alpha);
    else
        fprintf(fptr, "E %d %.5lf %.5lf\n", j+1, bestC[j].alpha, bestC[j].beta);
  }
  double threshold = bestC[n-1].alpha - 9.91152; // lhuang -9.xxx or ?
  for (j = 0; j < n; j++) 
    print_states(fptr, pfscore[j], j, "P", inside_only, threshold);
  for (j = 0; j < n; j++) 
    print_states(fptr, bestM[j], j, "M", inside_only, threshold);
  for (j = 0; j < n; j++) 
    print_states(fptr, bestM2[j], j, "M2", inside_only, threshold);
  for (j = 0; j < n; j++) 
    print_states(fptr, bestMulti[j], j, "Multi", inside_only, threshold);
}

BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose,
                             bool pf_allocate)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose),
      allocate(pf_allocate) {
      initialize();
}

BeamCKYParser::~BeamCKYParser(){
    delete[] nucs;
    delete[] bestC;
    delete[] bestH;
    // delete[] bestP;
    delete[] bestM;
    delete[] bestM2;
    delete[] bestMulti;
}