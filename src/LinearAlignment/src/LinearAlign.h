/*
 *LinearAlign.h*
 header file for LinearAlign.cpp.

 author: Sizhen Li
 created by: 06/2020
*/

#ifndef BEAM_ALIGNMENT_H
#define BEAM_ALIGNMENT_H

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <math.h> 
#include "../../utils/math/xlog_math.h"
#include "../../utils/math/matrix.h"
#include "../../utils/phmm.h"

using namespace std;

#define VALUE_MIN numeric_limits<double>::lowest()
#define GET_ACGU_NUM(x) ((x=='A'? 0 : (x=='C'? 1 : (x=='G'? 2 : (x=='U'?3: 4)))))
#define GET_NUC(x) ((x==0? 'A' : (x==1 ? 'C' : (x==2 ? 'G' : (x==3 ? 'U': '.')))))

struct match_score {
    double upstream;
	double downstream;
	double unpair;
	match_score(): upstream(0), downstream(0), unpair(0){};
};

struct AlignState {
	unsigned i;
	unsigned k;
    double alpha;
    double beta;
    int manner;
	int pre;
	unsigned step;
	bool close; 

    AlignState(): alpha(xlog(0)), beta(xlog(0)), manner(0), pre(0), close(false), step(0), i(0), k(0) {};

    void set(double score_, int pre_manner_, int manner_, unsigned step_, unsigned i_, unsigned k_) {
        alpha = score_; pre = pre_manner_; manner = manner_; step = step_ ; i = i_; k = k_;
    }

    void set(double score_, int manner_) {
        alpha = score_; manner = manner_;
    }

	void set_beta(double score_, int pre_manner_, int manner_, unsigned step_, unsigned i_, unsigned k_) {
        beta = score_; pre = pre_manner_; manner = manner_; step = step_ ; i = i_; k = k_;
    }
};

struct aln_ret {
    double prob = xlog(0.0);
    double aln_prob = xlog(0.0);
};

struct multi_aln_ret {
    float value = 0.0;
};

class LinearTurboFold;
class BeamAlign{
public:
    int beam;

    BeamAlign(int beam_size=100);
	~BeamAlign();
    LinearTurboFold* out;
    int i_seq1, i_seq2;

    void ml_alignment(string &seq1, string &seq2, vector<char> &aln1, vector<char> &aln2, double** &transprobs, double** &emitprobs, bool prior);
	double forward(string seq1, string seq2, double** &trans_probs, double** &emit_probs, bool prior);
	double backward(double** &transprobs, double** &emitprobs, bool prior);
	std::unordered_map<int, aln_ret>* cal_align_prob(double forward_score, double threshold, std::unordered_map<int, aln_ret>* &aln_ret);

private:
    unordered_map<int, AlignState> *bestINS1, *bestINS2, *bestALN; //, *bestState;
    int *nucs1, *nucs2;
    vector<pair<double, int>> scores;
    unsigned seq1_len, seq2_len, max_len;

    void prepare(string &len1, string &len2);
    double beam_prune(std::unordered_map<int, AlignState> &beamstep);
    double beam_prune(std::unordered_map<int, AlignState> &beamstep, unsigned i);
    double get_trans_emit_prob(int prev_state, int current_state, int i, int k, double** &transprobs, double** &emitprobs);
    void traceback(vector<char> &aln1, vector<char> &aln2);
    double get_match_prior(int i, int k, bool prior);
    double quickselect(vector<pair<double, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k);
    unsigned long quickselect_partition(vector<pair<double, int>>& scores, unsigned long lower, unsigned long upper);
};

#endif