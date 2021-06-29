/*
 *LinearPartition.h*
 header file for LinearPartition.cpp.

 author: He Zhang
 created by: 03/2019
*/

#ifndef FASTCKY_BEAMCKYPAR_H
#define FASTCKY_BEAMCKYPAR_H

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <math.h> 

#include "../../utils/pfunction_math.h"
#include "../../utils/structure.h"

// #define MIN_CUBE_PRUNING_SIZE 20
#define kT 61.63207755

#define NEG_INF -2e20 
// #define testtime

using namespace std;

#ifdef FAST_FLOAT
  typedef float pf_type;
#else
  typedef double pf_type;
#endif


#ifdef lpv
  typedef int value_type;
  #define VALUE_MIN numeric_limits<double>::lowest()
#else
  typedef float value_type;
  #define VALUE_MIN numeric_limits<double>::lowest()
#endif

enum Manner {
  MANNER_NONE = 0,              // 0: empty
  MANNER_H,                     // 1: hairpin candidate
  MANNER_HAIRPIN,               // 2: hairpin
  MANNER_SINGLE,                // 3: single
  MANNER_HELIX,                 // 4: helix
  MANNER_MULTI,                 // 5: multi = ..M2. [30 restriction on the left and jump on the right]
  MANNER_MULTI_eq_MULTI_plus_U, // 6: multi = multi + U
  MANNER_P_eq_MULTI,            // 7: P = (multi)
  MANNER_M2_eq_M_plus_P,        // 8: M2 = M + P
  MANNER_M_eq_M2,               // 9: M = M2
  MANNER_M_eq_M_plus_U,         // 10: M = M + U
  MANNER_M_eq_P,                // 11: M = P
  MANNER_C_eq_C_plus_U,         // 12: C = C + U
  MANNER_C_eq_C_plus_P,         // 13: C = C + P
};

// A hash function used to hash a pair of any kind 
struct hash_pair { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1, T2>& p) const
    { 
        auto hash1 = hash<T1>{}(p.first); 
        auto hash2 = hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
};

struct ExtValue {
    double score;
    ExtValue(): score(LOG_OF_ONE) {};
};

struct State {
    float alpha;
    float beta;
    float score_external;

    State(): alpha(VALUE_MIN), beta(VALUE_MIN) {};
};

class LinearTurboFold;
class BeamCKYParser {
public:
    int beam;
    bool no_sharp_turn;
    bool is_verbose;
    string bpp_file;
    string forest_file;
   
    struct DecoderResult {
        float alpha;
        // unsigned long num_states;
        double time;
    };

    BeamCKYParser(int beam_size=100,
                  bool nosharpturn=true,
                  bool is_verbose=false,
                  bool allocate=false);
    ~BeamCKYParser();

    double parse(int i_iter, int i_seq, LinearTurboFold* outer, string& seq, unordered_map<int, State>* &pfscore, unordered_map<int, ExtValue>* &ext_info, string pfSaveFile, string bppSaveFile);
    void getbasepairs(int j, vector<int>& pairs);
    bool allocate;

private:
    void get_parentheses(char* result, string& seq);

    unsigned seq_length;

    unordered_map<int, State> *bestH, *bestM2, *bestMulti, *bestM; // *bestP -> pfscore
 
    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    State *bestC;

    int *nucs;

    void prepare(unsigned len);


    void cal_PairProb(State& viterbi, unordered_map<int, State>* &pfscore, unordered_map<int, ExtValue>* &ext_info); 

    void outside(LinearTurboFold* outer, int i_seq, int i_iter, vector<int> next_pair[], unordered_map<int, State>* &pfscore, unordered_map<int, ExtValue>* &ext_info);

    float beam_prune(unordered_map<int, State>& beamstep);
    float beam_prune(int i_seq, LinearTurboFold* outer, std::unordered_map<int, State> &beamstep, int j, unordered_map<int, ExtValue>* &ext_info);

    vector<pair<float, int>> scores;

    unordered_map<pair<int,int>, float, hash_pair> Pij;

    void output_to_file(string file_name, const char * type);
    void print_states(FILE *fptr, unordered_map<int, State>& states, int j, string label, bool inside_only, double threshold);
    void dump_forest(string seq, unordered_map<int, State>* &pfscore, bool inside_only);

};


// log space: borrowed from CONTRAfold

inline float Fast_LogExpPlusOne(float x){
  
    // Bounds for tolerance of 7.05e-06: (0, 11.8625)
    // Approximating interval: (0, 0.661537) --> ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306);
    // Approximating interval: (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989);
    // Approximating interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    // Approximating interval: (2.49126, 3.37925) --> ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829);
    // Approximating interval: (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399);
    // Approximating interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    // Approximating interval: (5.78907, 7.81627) --> ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903);
    // Approximating interval: (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051);
    // 8 polynomials needed.
    
    assert(float(0.0000000000) <= x && x <= float(11.8624794162) && "Argument out-of-range.");
    if (x < float(3.3792499610))
    {
        if (x < float(1.6320158198))
        {
            if (x < float(0.6615367791))
                return ((float(-0.0065591595)*x+float(0.1276442762))*x+float(0.4996554598))*x+float(0.6931542306);
            return ((float(-0.0155157557)*x+float(0.1446775699))*x+float(0.4882939746))*x+float(0.6958092989);
        }
        if (x < float(2.4912588184))
            return ((float(-0.0128909247)*x+float(0.1301028251))*x+float(0.5150398748))*x+float(0.6795585882);
        return ((float(-0.0072142647)*x+float(0.0877540853))*x+float(0.6208708362))*x+float(0.5909675829);
    }
    if (x < float(5.7890710412))
    {
        if (x < float(4.4261691294))
            return ((float(-0.0031455354)*x+float(0.0467229449))*x+float(0.7592532310))*x+float(0.4348794399);
        return ((float(-0.0010110698)*x+float(0.0185943421))*x+float(0.8831730747))*x+float(0.2523695427);
    }
    if (x < float(7.8162726752))
        return ((float(-0.0001962780)*x+float(0.0046084408))*x+float(0.9634431978))*x+float(0.0983148903);
    return ((float(-0.0000113994)*x+float(0.0003734731))*x+float(0.9959107193))*x+float(0.0149855051);

    /*
    // Bounds for tolerance of 9.99e-05: (0, 9.21129)
    // Approximating interval: (0, 1.40131) --> ((T(-0.0118287252)*x+T(0.1342168806))*x+T(0.4976005362))*x+T(0.6932470806);
    // Approximating interval: (1.40131, 3.06792) --> ((T(-0.0117040733)*x+T(0.1232945547))*x+T(0.5276092444))*x+T(0.6721240615);
    // Approximating interval: (3.06792, 5.15409) --> ((T(-0.0027005983)*x+T(0.0419040665))*x+T(0.7762991688))*x+T(0.4152395732);
    // Approximating interval: (5.15409, 9.21129) --> ((T(-0.0001617326)*x+T(0.0040111354))*x+T(0.9666890441))*x+T(0.0929363811);
    // 4 polynomials needed.
    
    assert(float(0.0000000000) <= x && x <= float(9.2112909219), "Argument out-of-range.");
    if (x < float(3.0679202382))
    {
        if (x < float(1.4013117629))
            return ((float(-0.0118287252)*x+float(0.1342168806))*x+float(0.4976005362))*x+float(0.6932470806);
        return ((float(-0.0117040733)*x+float(0.1232945547))*x+float(0.5276092444))*x+float(0.6721240615);
    }
    if (x < float(5.1540922927))
        return ((float(-0.0027005983)*x+float(0.0419040665))*x+float(0.7762991688))*x+float(0.4152395732);
    return ((float(-0.0001617326)*x+float(0.0040111354))*x+float(0.9666890441))*x+float(0.0929363811);
    */
}

inline void Fast_LogPlusEquals (float &x, float y)
{
    if (x < y) std::swap (x, y);
    if (y > float(NEG_INF/2) && x-y < float(11.8624794162))
        x = Fast_LogExpPlusOne(x-y) + y;
}

inline float Fast_Exp(float x)
{
    // Bounds for tolerance of 4.96e-05: (-9.91152, 0)
    // Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
    // Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
    // Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
    // Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
    // Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
    // Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
    // 6 polynomials needed.
    
    if (x < float(-2.4915033807))
    {
        if (x < float(-5.8622823336))
        {
            if (x < float(-9.91152))
                return float(0);
            return ((float(0.0000803850)*x+float(0.0021627428))*x+float(0.0194708555))*x+float(0.0588080014);
        }
        if (x < float(-3.8396630909))
            return ((float(0.0013889414)*x+float(0.0244676474))*x+float(0.1471290604))*x+float(0.3042757740);
        return ((float(0.0072335607)*x+float(0.0906002677))*x+float(0.3983111356))*x+float(0.6245959221);
    }
    if (x < float(-0.6725053211))
    {
        if (x < float(-1.4805375919))
            return ((float(0.0232410351)*x+float(0.2085645908))*x+float(0.6906367911))*x+float(0.8682322329);
        return ((float(0.0573782771)*x+float(0.3580258429))*x+float(0.9121133217))*x+float(0.9793091728);
    }
    if (x < float(0))
        return ((float(0.1199175927)*x+float(0.4815668234))*x+float(0.9975991939))*x+float(0.9999505077);
    return (x > float(46.052) ? float(1e20) : expf(x));
}

#endif //FASTCKY_BEAMCKYPAR_H
