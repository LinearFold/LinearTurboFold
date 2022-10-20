/*
 *bpp.h*
 The main code for base pair probability calculation.

 author: He Zhang
 created by: 04/2019
*/

#include <stdio.h> 
#include <sys/time.h>
#include "LinearPartition.h"
#include "Utils/utility.h"
#include "../../LinearTurboFold.h"
// #define lpv

using namespace std;

void BeamCKYParser::output_to_file(string file_name, const char * type) {

    if(!file_name.empty()) {
        // printf("Outputing base pairing probability matrix to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }

        int turn = no_sharp_turn?3:0;
        for (int i = 1; i <= seq_length; i++) {
            for (int j = i + turn + 1; j <= seq_length; j++) {
                pair<int, int> key = make_pair(i,j);
                auto got = Pij.find(key);
                if (got != Pij.end()){
                    fprintf(fptr, "i=%d, j=%d, probs=%.4e\n", i, j, got->second);
                }
            }
        }
        fprintf(fptr, "\n");
        fclose(fptr); 
        // printf("Done!\n"); 
    }

    return;
}

void BeamCKYParser::cal_PairProb(State& viterbi, unordered_map<int, State>* &pfscore, unordered_map<int, ExtValue>* &ext_info) {

    for(int j=0; j<seq_length; j++){
        for(auto &item : pfscore[j]){
            int i = item.first;
            State state = item.second;
            
            float temp_prob_inside = state.alpha + state.beta - viterbi.alpha - ext_info[j+1][i+1].score;
            if (temp_prob_inside > float(-9.91152)) {
                float prob = Fast_Exp(temp_prob_inside);
                if(prob > 1.0) prob = 1.0;
                // if(prob < bpp_cutoff) continue;
                Pij[make_pair(i+1, j+1)] = prob;
            }
        }
    }

    // output to a single file with user specified name;
    if (!bpp_file.empty()){
        output_to_file(bpp_file, "w");
    }

    return;
}

void BeamCKYParser::outside(LinearTurboFold* outer, int i_seq, int i_iter, vector<int> next_pair[], unordered_map<int, State>* &pfscore, unordered_map<int, ExtValue>* &ext_info){
      
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    bestC[seq_length-1].beta = 0.0;

    // from right to left
    value_type newscore;
    for(int j = seq_length-1; j >= 0; --j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = pfscore[j]; // bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
                Fast_LogPlusEquals(beamstepC.beta, (bestC[j+1].beta));
            }
        }

        if (j == 0) continue;
    
        // beam of M
        {
            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta);
                }
            }
        }

        // beam of M2
        {
            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];
                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
                            Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta);
                        }
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(state.beta, beamstepM[i].beta);
            }
        }

        // beam of P
        {  
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                if (i >0 && j<seq_length-1) {
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            // lisiz: if (bestP[q].find(p) != bestP[q].end()) {
                            if (pfscore[q].find(p) != pfscore[q].end()) {
                                if (p == i - 1 && q == j + 1) {
                                    // helix
                                    int score_single = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                                nuci_1, nuci, nucj, nucj1);
                                    Fast_LogPlusEquals(state.beta, (pfscore[q][p].beta + score_single/kT));
                                    // Fast_LogPlusEquals(state.beta, (bestP[q][p].beta + score_single/kT));
                                } else {
                                    // single branch
                                    int score_single = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                    nuci_1, nuci, nucj, nucj1);
                                    Fast_LogPlusEquals(state.beta, (pfscore[q][p].beta + score_single/kT));
                                    // Fast_LogPlusEquals(state.beta, (bestP[q][p].beta + score_single/kT));
                                }
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }

                // 2. M = P
                if(i > 0 && j < seq_length-1){
                        int score_M1 = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, (beamstepM[i].beta + score_M1/kT));
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
                    int M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    float m1_alpha = M1_score/kT;
                    float m1_plus_P_alpha = state.alpha + m1_alpha;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(state.beta, (beamstepM2[newi].beta + m_state.alpha + m1_alpha));
                        Fast_LogPlusEquals(m_state.beta, (beamstepM2[newi].beta + m1_plus_P_alpha));
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        int nuck = nuci_1;
                        int nuck1 = nuci;
                        int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                 nucj, nucj1, seq_length);
                        float external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + score_external_paired/kT;
                        Fast_LogPlusEquals(bestC[k].beta, state.alpha + external_paired_alpha_plus_beamstepC_beta);
                        Fast_LogPlusEquals(state.beta, bestC[k].alpha + external_paired_alpha_plus_beamstepC_beta);
                    } else {
                        // value_type newscore;
                        int score_external_paired = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, (beamstepC.beta + score_external_paired/kT));
                    }
                }

                if (i_iter > 0) state.beta = PROD(state.beta, ext_info[j+1][i+1].score);
            }
        }

        // beam of Multi
        {
            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
                        Fast_LogPlusEquals(state.beta, (bestMulti[jnext][i].beta));
                    }
                }

                // 2. generate P (i, j)
                {
                    int score_multi = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    if (beamstepP.find(i) == beamstepP.end()) continue;
                    Fast_LogPlusEquals(state.beta, (beamstepP[i].beta + score_multi/kT));
                }
            }
        }
    }  // end of for-loo j

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    return;
}

