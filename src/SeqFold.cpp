#include "SeqFold.h"
#include "probknot.h"

// Copy Thermodynamics info and datatable from another Thermodynamics instance.
SeqFold::SeqFold(int ckyBeam, const char * filepathOrSequence, bool verbose) {
    // allocate ct
    ct = new structure();

    // Indicate that the partition function calculation has not been performed.
    partitionfunctionallocated = false;

    // Load the sequence or file
    if (filepathOrSequence == NULL) return; 
    ct->SetSequence(filepathOrSequence);

    // LinearPartition
    bool sharpturn = false; // TODO
    linearpf = new BeamCKYParser(ckyBeam, !sharpturn, verbose);
}

// lisiz
// Calculate the partition function using LinearPartition
int SeqFold::LinearPartition(int i_iter, int i_seq, string pfsSaveFile, string bppSaveFile){
    int i, j;

    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;
    
    //The thermodynamic data tables have not been read 
    //if (!VerifyThermodynamic()) return 5; 

    // Allocate the memory needed (only if this is the first call to pfunction):
    // indicate that the memory has been allocated so that the destructor will delete it.
    partitionfunctionallocated = true;

    string seq = string(ct->GetSequence());
    
    // create
    pfscore = new unordered_map<int, State>[ct->GetSequenceLength()];
    extrisic_info = new unordered_map<int, ExtValue>[ct->GetSequenceLength() + 1];
    viterbi_score = linearpf->parse(i_iter, i_seq, this->outer, seq, pfscore, extrisic_info, pfsSaveFile, bppSaveFile);

    return 0;
}

// lisiz
void SeqFold::Savepartitionfunction(int i_iter){
    // TODO
    int i, j;

    if (i_iter > 0) {
        delete[] prepfscore;
        delete[] pre_ext_info;
    } 
    // save bestP
    prepfscore = new unordered_map<int, State>[ct->GetSequenceLength()];
    for(int j=0; j<ct->GetSequenceLength(); j++){
        // prepfscore[j] = std::move(pfscore[j]); //
        // cout << j << " size: " << pfscore[j].size() << endl;
        for(auto &item : pfscore[j]){
            int i = item.first;
            State state = item.second;
            
            prepfscore[j][i].alpha = state.alpha;
            prepfscore[j][i].beta = state.beta;
        }
    }
    delete[] pfscore;
    previterbi = viterbi_score;

    // save extrinsic info
    pre_ext_info = new unordered_map<int, ExtValue>[ct->GetSequenceLength() + 1];
    for(int j=1; j<=ct->GetSequenceLength(); j++){
        for(auto &item : extrisic_info[j]){
            int i = item.first;
            pre_ext_info[j][i].score = item.second.score;
        }
    }
    delete[] extrisic_info;

    // char *savefilename;
    // if (is_blank(savefile)) 
	//     saveFile=NULL;
    // else {
    //     saveFile=new char[((int) strlen(savefile))+1];
    //     strcpy(saveFile, savefile);
    // }
}

void SeqFold::ComputeMatchScore(match_score** &match_score_ret, int i_seq){
    for (int i = 0; i < ct->GetSequenceLength(); i++){  // prepfscore start from 0
        for (auto &state : prepfscore[i]){
            int k = state.first;
            double score = GetBasePairProb(k+1, i+1);
            match_score_ret[i_seq][i].upstream += score;
            match_score_ret[i_seq][k].downstream += score;
            // if ((i_seq == 0) and ((i == 2) or (k == 2))){
            //     printf("%d, %d, %f\n", i, k , score);
            // }
        }
    }
    for (int i = 0; i < ct->GetSequenceLength(); i++){
        double upscore = match_score_ret[i_seq][i].upstream;
        double downscore = match_score_ret[i_seq][i].downstream;
        double unpair = 1 - upscore - downscore;
        // if (i < 3){
        //     printf("%d, %d, %f, %f, %f\n", i_seq, i, upscore, downscore, unpair);
        // }
        if (unpair < 0.0){
            if (upscore > 1.){
                upscore = 1.0;
                downscore = 0.0;
            } else if (downscore > 1.)
            {
                upscore = 0.0;
                downscore = 1.0;
            }
            unpair = 0.0;
        }
        match_score_ret[i_seq][i].unpair = unpair;
        match_score_ret[i_seq][i].upstream = upscore;
        match_score_ret[i_seq][i].downstream = downscore;
        // cout << i_seq << " " << i << " " << upscore << " " << downscore << " " << unpair << endl;
    }
}

double SeqFold::GetBasePairProb(int i, int j){
    if (i==0 || j==0) return 0.0;
    if (prepfscore[j-1].find(i-1) == prepfscore[j-1].end()) return 0.0; // N.B. otherwise increase memory usage
    if (pre_ext_info[j][i].score < TO_XLOG(EPSILON)) return 0.0; // -690.776

    // double prob = DIV(PROD(v->f(i,j), v->f(j,i+ct->GetSequenceLength())), PROD(w5[ct->GetSequenceLength()] POWSCALING2(2), ct->extrinsic_info[j][i]));
    // if ((i == 3) or (j == 3)){
    //     printf("%f, %f, %f, %f\n", prepfscore[j-1][i-1].alpha, prepfscore[j-1][i-1].beta, previterbi, pre_ext_info[j][i].score);
    // }
    double prob = DIV(PROD(prepfscore[j-1][i-1].alpha, prepfscore[j-1][i-1].beta), PROD(previterbi POWSCALING2(2), pre_ext_info[j][i].score));
    if (prob > float(-9.91152)) {
        // if (prob > 0.1) cout << i << " " << j << " " << v->f(i,j) << " " <<  v->f(j,i+ct->GetSequenceLength()) << " " << w5[ct->GetSequenceLength()]<< " " << ct->extrinsic_info[j][i] << " " << prob << endl;
        prob = TO_LINEAR(prob);
        if (prob > 1.0) prob = 1.0; // lisiz: no pairs with prob > 1.0, TBD
        return prob;
    }
    return 0.0;
}

structure *SeqFold::GetStructure() {
    return ct;
}

void SeqFold::SetSequenceLabel(const string& label) {
    GetStructure()->SetSequenceLabel(label);
}

int SeqFold::GetErrorCode() const {
    return ErrorCode;
}

SeqFold::~SeqFold() {
    if (partitionfunctionallocated) {
        // The partition function calculation was performed, so the memory allocated for the partition function needs to be deleted.
        delete[] prepfscore;
        delete[] pre_ext_info;
        delete[] saveFile;
    }

    delete ct;//delete the structure

    delete linearpf;
}

//This is a protected function for handling file input.
int SeqFold::FileReader(const char filename[], const RNAInputType type) {

    if (!isStdIoFile(filename) && !fileExists(filename)) {
        // SetErrorDetails(sfmt("The path '%s' is invalid or does not exist.", filename));
        cout << "The path '%s' is invalid or does not exist." << endl;
        return 1; // file not found.
    }

    // if (type==FILE_CT||type==FILE_SEQ||type==FILE_DBN)
    //     if (!IsAlphabetRead()) return 30;

    // RMW 2015-03-12: try/catch to prevent any unexpected errors from crashing the program. At least show an error message.
    // (Previous to this, passing the wrong type of file to openct caused an unhandled memory allocation exception.)
    try {
        //open the file based on type:
        switch(type) {
            // case FILE_CT: // ct file
            //     return ct->openct(filename);  // openct returns 0 on success and uses the same non-zero error codes as RNA::GetErrorMessage. So the result can be returned directly.
            // case FILE_DBN: // dot-bracket file
            //     return ct->opendbn(filename); // opendbn returns 0 on success and uses the same non-zero error codes as RNA::GetErrorMessage. So the result can be returned directly.
            case FILE_SEQ:
                //type indicates a .seq file
                return ct->openseqx(filename);
            default:
                return 22; // error - invalid file type
        } // SWITCH type
    } catch (std::exception* ex) {
        // SetErrorDetails(ex->what());
        cout << "can not read file" << endl;
        return 2;
    }
}

int SeqFold::ProbKnot(int iterations, int MinHelixLength, double threshold) {
    if (iterations < 1) {
        printf("there can't be fewer than one iteration"); 
        return 0;
    }

	if (threshold < 0) {
		printf("There can't be pairs with less than zero probablity");
		return 0; 
	}

    //Past error trapping
    //Call the ProbKnot Program:
    return LinearProbKnotAssemble(ct, prepfscore, pre_ext_info, previterbi, iterations, MinHelixLength, threshold);
}

int SeqFold::WriteCt(const char filename[], bool append, CTCommentProvider &commentProvider) const {
    if (ct->GetNumberofStructures()>0)
        return ct->ctout(filename,append,commentProvider);
    
    return 0;
}

//! If there was an error, this returns the error message, along with error details (if any).
//! if there was no error (i.e. GetErrorCode() returns 0 and GetErrorDetails() returns NULL), this function returns an empty string ("");
// string SeqFold::GetFullErrorMessage() const {
//     int code = GetErrorCode();
//     string message(code==0?"":GetErrorMessage(code));
//     string details = GetErrorDetails();

//     // If message and details are both non-empty, 
//     // combine them to "<message>: <details>"
//     // This requires triming whitepace and the dot/period (.) from the end of message.
//     if (!message.empty() && !details.empty()) {
//         std::size_t last = message.find_last_not_of("\r\n\t .");
//         if (last != string::npos)
//             message.resize(last+1);
//         message.append(": ");
//     }
//     message.append(details);
//     // GetErrorMessage always ends with \n, so this should also, for consistency.
//     if (!message.empty() && message[message.length()-1]!='\n')
//         message+='\n';
//     return message;
// }
