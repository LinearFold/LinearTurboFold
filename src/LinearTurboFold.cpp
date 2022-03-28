/*
 *LinearSampling.cpp*
 The main code for LinearTurboFold: Linear-Time RNA Structural
 Alignment and Conserved Structure Prediction
 with Applications in Coronaviruses.

 author: Sizhen Li
 edited by: 10/2020
*/

#include <fstream>
#include <iostream>
#include <unistd.h>
#include <sys/time.h>
#include "LinearTurboFold.h"

using namespace std;

double LinearTurboFold::get_folding_extrinsic_information(int i_seq1, int j, int i){
    double seq1_seq2_seq_similarity_weight;
    double seq1_seq2_mapping_probability;
    double extrinsic_info = 0.0f;

    double norm = 0.0;
    for (unsigned int i_seq2 = 0; i_seq2 < sequences.size(); i_seq2++){
        if (i_seq2 == i_seq1) continue;

        if (i_seq1 < i_seq2) seq1_seq2_seq_similarity_weight = (1.0f - similarities[i_seq1][i_seq2]);
        else seq1_seq2_seq_similarity_weight = (1.0f - similarities[i_seq2][i_seq1]);

        norm += 2 * seq1_seq2_seq_similarity_weight;
    }

    for (unsigned int i_seq2 = 0; i_seq2 < sequences.size(); i_seq2++){
        if (i_seq2 == i_seq1) continue;
        SeqFold* RNA2 = folds[i_seq2];

        if (i_seq1 < i_seq2) seq1_seq2_seq_similarity_weight = (1.0f - similarities[i_seq1][i_seq2]);
        else seq1_seq2_seq_similarity_weight = (1.0f - similarities[i_seq2][i_seq1]);

        for(auto &i_item : beam_aln_results[i_seq1][i_seq2][i]){
            int k = i_item.first;
            double ikprob = i_item.second.prob;

            for(auto &j_item : beam_aln_results[i_seq1][i_seq2][j]){
                int l = j_item.first;
                if (k == 0 || l == 0) continue;
                double jlprob = j_item.second.prob;
            
                seq1_seq2_mapping_probability = ikprob * jlprob;
                extrinsic_info += 2 * seq1_seq2_seq_similarity_weight * seq1_seq2_mapping_probability * RNA2->GetBasePairProb(k, l); // why two times??
            }
        } 
    }

    return TO_XLOG(extrinsic_info); // / norm);
}

LinearTurboFold::LinearTurboFold(vector<t_structure*> *fasta_sequences, int ckyBeam, int hmmBeam, bool verbose){
    err_code = 0;
    progress=NULL;
    multiple_sequences=NULL;
    multiple_alignment=NULL;

    hmm_beam = hmmBeam;
    cky_beam = ckyBeam;
    isVerbose = verbose;

    num_seq = fasta_sequences->size();
    if(num_seq == 0) {
        cout << "Need at least 1 sequence to predict structure for." << endl; // TODO: error output
        return;
    }
    sequences.resize(num_seq);
    folds.resize(num_seq);
    saves.resize(num_seq);

    for(unsigned int i_str = 0; i_str < fasta_sequences->size(); i_str++)
    {
        sequences[i_str] = (*fasta_sequences)[i_str]; //new t_structure(fasta_sequences[i_str]);
        SeqFold* seqfold = new SeqFold(ckyBeam, &(sequences[i_str]->nucs[1]), verbose);
        seqfold->outer = this; // call get ext info function in LP
        seqfold->SetSequenceLabel(sequences[i_str]->ctlabel);
        // Now that CT label has been read, we can replace characters that are invalid in alignments.
        sequences[i_str]->check_set_label();
        folds[i_str]=seqfold;
    } 

    // initialize alignment information
    allocate_phmm();
    initialize_multiple_sequences();
}

void LinearTurboFold::run_iterations(ConfigParser config){
    n_iterations = config.turboIterations;
    struct timeval parse_starttime, parse_endtime;

    // Initialize the loops.
    for(int i_iter = 0; i_iter <= n_iterations; i_iter++)
    {
        // Add a coarse update of progress:
        if (progress!=NULL) {
            progress->update((int)((100.0*((double) i_iter))/((double) n_iterations+1)));
        }
        

        if (isVerbose){
            gettimeofday(&parse_starttime, NULL);
            printf("\n");
        }
        
        // pairwise alignments
        if (i_iter > 0) run_phmm_alignment(i_iter); 
        // each sequence folding
        refoldSequences(i_iter, config.saveBpps, config.savePfs, &config.outputPfsFiles, &config.outputBppFiles);

        if (isVerbose){
            gettimeofday(&parse_endtime, NULL);
            double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
            printf("Iter %d Done (%f seconds)\n", i_iter, parse_elapsed_time);
        }

        // Output final alignment
        if(i_iter == n_iterations)
        {   
            if (isVerbose) gettimeofday(&parse_starttime, NULL);
            run_phmm_alignment(i_iter);
            run_multiple_alignment();
            if (isVerbose) {
                gettimeofday(&parse_endtime, NULL);
                double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
                printf("MSA Generation Done (%f seconds)\n", parse_elapsed_time);
            }

        }
    } // i_iter loop.

    // Add a coarse update of progress:
    if (progress!=NULL) {
        progress->update(100);
    }

    // Printing info if write PFs or BPPs to files
    printf("%d iterations Done!\n", n_iterations);
    if (config.savePfs)
        printf("Outputing partition functions to files ...\n");
    if (config.saveBpps)
        printf("Outputing base pair probabilities to files ...\n");
    
    // Write Alignment file
    if (!config.OutAln.empty()) {
        printf("Outputing multiple sequence alignment to %s...\n", config.OutAln.c_str()); 
        FILE *fptr = fopen(config.OutAln.c_str(), "w"); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }

        bool isClustal = config.AlnFormat != "Fasta";
        multiple_alignment->WriteALN (fptr, config.ColumnNumber, isClustal);

        fprintf(fptr, "\n");
        fclose(fptr); 
    }

    // Secondary structure prediction and save structures
    printf("Outputing structures to files ...\n");
    for(int i_seq = 1; i_seq <= num_seq; i_seq++) {
        int ret = folds[i_seq-1]->ProbKnot(config.pkIterations, config.minHelixLength, config.threshold);  // TODO: parameters
        // printf("Outputing structures to %s\n", config.outputCtFiles[i_seq-1].c_str());
        int writeError = WriteCt(i_seq-1, config.outputCtFiles[i_seq-1].c_str());

        // ct2dot
        folds[i_seq-1]->GetStructure()->writedotbracket(config.outputDotFiles[i_seq-1].c_str());
    }

    return;
}


int LinearTurboFold::run_multiple_alignment()
{
    SafeVector<SafeVector<float> > distances (this -> GetNumberSequences(), SafeVector<float> (this -> GetNumberSequences(), 0));
    ProbabilisticModel model;
    vector<vector<unordered_map<int, multi_aln_ret>*>> mul_aln_results;
    mul_aln_results.resize(sequences.size());
    for(unsigned int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++){
        mul_aln_results[i_seq1].resize(sequences.size());
        for(unsigned int i_seq2 = 0; i_seq2 < sequences.size(); i_seq2++) {
            if(i_seq1 == i_seq2) continue;
            mul_aln_results[i_seq1][i_seq2] = new unordered_map<int, multi_aln_ret>[sequences[i_seq1]->numofbases + 1];
        }
    }
    for(unsigned int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
    {
        for(unsigned int i_seq2 = i_seq1+1; i_seq2 < sequences.size(); i_seq2++)
        {
            if(i_seq1 != i_seq2)
            {   
                for(int i = 0; i <= sequences[i_seq1]->numofbases; i++) {
                    for(auto &item : beam_aln_results[i_seq1][i_seq2][i]){
                        int k = item.first;
                        double value = item.second.aln_prob;
                        if (value >= 0.01) {
                            mul_aln_results[i_seq1][i_seq2][i][k].value = beam_aln_results[i_seq1][i_seq2][i][k].aln_prob;
                            mul_aln_results[i_seq2][i_seq1][k][i].value = beam_aln_results[i_seq1][i_seq2][i][k].aln_prob;
                        }
                    }
                }

                const unordered_map<int, multi_aln_ret>* const_mul_aln_results = mul_aln_results[i_seq1][i_seq2];
                pair<SafeVector<char> *, float> pair_alignment = model.LinearComputeAlignment(hmm_beam, sequences[i_seq1]->numofbases,sequences[i_seq2]->numofbases, const_mul_aln_results);

                float distance = pair_alignment.second / min (sequences[i_seq1]->numofbases,sequences[i_seq2]->numofbases);
                distances[i_seq1][i_seq2] = distances[i_seq2][i_seq1] = distance;
                delete pair_alignment.first;
            }
        } // i_seq2 loop.
    } // i_seq1 loop.


    for (int r = 0; r<numConsistencyReps; r++ ) {
        vector<vector<unordered_map<int, multi_aln_ret>*>> new_aln_results = LinearMultiConsistencyTransform(this->multiple_sequences, mul_aln_results);
        // now replace the old posterior probs.
        for (int i = 0; i < this -> GetNumberSequences(); i++){
            for (int j = 0; j < this -> GetNumberSequences(); j++){
                if (i == j) continue;
                delete[] mul_aln_results[i][j];
                mul_aln_results[i][j] = new_aln_results[i][j];
            }
        }
    }
    
    this->multiple_sequences->SaveOrdering();

    this->multiple_alignment=NULL;
    TreeNode *tree = TreeNode::ComputeTree(distances); // lisiz, guide tree 
  
    // make the final alignment
    this->multiple_alignment = LinearComputeFinalAlignment(tree, this->multiple_sequences, mul_aln_results, model, hmm_beam);
    int numSeqs = this->multiple_sequences->GetNumSequences();
    for (int i = 0; i < numSeqs; i++){
        mul_aln_results[i].clear();
    }
    mul_aln_results.clear();

    delete tree;

    return(0);
}

// refold all sequences.
void LinearTurboFold::refoldSequences(int i_iter, bool saveBpps, bool savePfs, vector<string>* pfSaveFiles, vector<string>* bppSaveFiles){
    for(int i_seq = 0; i_seq < num_seq; i_seq++) {
        SeqFold &seqfold = *folds[i_seq];

        string pfSaveFile = "";
        string bppSaveFile = "";
        // On the last iteration, save PFs and BPPs if the user has specified one.
        if ((i_iter == n_iterations) && (savePfs))
            pfSaveFile = (*pfSaveFiles)[i_seq];
        if ((i_iter == n_iterations) && (saveBpps))
            bppSaveFile = (*bppSaveFiles)[i_seq];
        
        seqfold.LinearPartition(i_iter, i_seq, pfSaveFile, bppSaveFile); 
    }

    // save partition function for the next iteration
    // calculate match score for the alignment computation in the next iteration
    for (unsigned int i_seq = 0; i_seq < sequences.size(); i_seq++){
        SeqFold &seqfold = *folds[i_seq];
        seqfold.Savepartitionfunction(i_iter);
        seqfold.ComputeMatchScore(match_score_rets, i_seq);
    }
}

void LinearTurboFold::initialize_multiple_sequences()
{
    this->multiple_sequences=new MultiSequence();
    for(int i = 0 ; i < this->GetNumberSequences(); i++){
        SafeVector<char>* data = new SafeVector<char>(1+sequences[i]->numofbases);
        (*data)[0]='@';
        for (int b=1; b<=sequences[i]->numofbases;b++){
            (*data)[b]=toupper(sequences[i]->nucs[b]);
        }
        Sequence* temp_seq=new Sequence(data,string(sequences[i]->ctlabel),sequences[i]->numofbases,i,i);
        this->multiple_sequences->AddSequence(temp_seq);
    }
}

int LinearTurboFold::allocate_phmm()
{
    similarities = (double**)malloc(sizeof(double*) * (sequences.size() + 2));
    beam_aln_results.resize(sequences.size());
    match_score_results.resize(sequences.size());
    match_score_rets = new match_score*[sequences.size()];

    for(unsigned int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
    {
        similarities[i_seq1] = (double*)malloc(sizeof(double) * (sequences.size() + 2));
        match_score_rets[i_seq1] = new match_score[sequences[i_seq1]->numofbases + 1];
        beam_aln_results[i_seq1].resize(sequences.size());

        for(unsigned int i_seq2 = i_seq1+1; i_seq2 < sequences.size(); i_seq2++)
        {
            similarities[i_seq1][i_seq2] = 0.0f;                    
            if(i_seq1 != i_seq2)
            {
                beam_aln_results[i_seq1][i_seq2] = new unordered_map<int, aln_ret>[sequences[i_seq1]->numofbases + 1];
            }
        }
    }

    return(0);
}

int LinearTurboFold::run_phmm_alignment(int i_iter)
{
    bool using_prior = false;
    if(i_iter > 0) using_prior = true;

    for(unsigned int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
    {
        for(unsigned int i_seq2 = i_seq1+1; i_seq2 < sequences.size(); i_seq2++)
        {
            if(i_seq1 != i_seq2)
            {
                t_phmm_aln *phmm_aln = new t_phmm_aln(sequences[i_seq1], sequences[i_seq2]);
	            linearAlign = new BeamAlign(hmm_beam);
                linearAlign->out = this;
                linearAlign->i_seq1 = i_seq1;
                linearAlign->i_seq2 = i_seq2;
                phmm_aln->phmm = new t_phmm(ML_emit_probs, ML_trans_probs);

                if(i_iter > 0)
                {
                    delete[] beam_aln_results[i_seq1][i_seq2];
                    delete[] beam_aln_results[i_seq2][i_seq1];
                    beam_aln_results[i_seq1][i_seq2] = new unordered_map<int, aln_ret>[sequences[i_seq1]->numofbases + 1];
                    beam_aln_results[i_seq2][i_seq1] = new unordered_map<int, aln_ret>[sequences[i_seq2]->numofbases + 1];
                } else {
                    beam_aln_results[i_seq2][i_seq1] = new unordered_map<int, aln_ret>[sequences[i_seq2]->numofbases + 1];
                }   

                string seq1(sequences[i_seq1]->linearnucs); // constructor
                string seq2(sequences[i_seq2]->linearnucs);
                vector<char> aln1, aln2;
                linearAlign->ml_alignment(seq1, seq2, aln1, aln2, phmm_aln->phmm->trans_probs, phmm_aln->phmm->emission_probs, using_prior);
                delete(phmm_aln->phmm); // This is not needed after the computation any more.

                int aln_len = aln1.size();
                assert (aln1.size() == aln2.size());

                // print alignment
                char* seq1_aln_line_str = (char*)malloc(sizeof(char) * (aln_len+1)); 
                char* seq2_aln_line_str = (char*)malloc(sizeof(char) * (aln_len+1));

                for (int i = aln_len - 1; i >= 0; --i){
                    seq1_aln_line_str[aln_len - 1 - i] = aln1[i];
                    seq2_aln_line_str[aln_len - 1 - i] = aln2[i];
                }
                seq1_aln_line_str[aln_len] = 0;
                seq2_aln_line_str[aln_len] = 0;
                // cout << seq1_aln_line_str << endl;
                // cout << seq2_aln_line_str << endl;

                t_p_alignment* ml_aln = new t_p_alignment(seq1_aln_line_str, seq2_aln_line_str);
                double similarity = ml_aln->get_aln_similarity(GAP_SYM);
                similarities[i_seq1][i_seq2] = similarity;
                
                free(seq1_aln_line_str);
                free(seq2_aln_line_str);
                delete(ml_aln);

                // Validate and load parameter files: The loop constants. 
                char *data_dir = getcwd(NULL, 0);
                if(data_dir == NULL)
                {
                    printf("Could not resolve thermodynamics data directory.\n");
                    exit(0);
                }
                char* phmm_pars_fp = (char*)malloc(sizeof(char*) * (strlen(data_dir) + strlen("src/data_tables/fam_hmm_pars.dat") + 2));
                sprintf(phmm_pars_fp, "%s/%s", data_dir, "src/data_tables/fam_hmm_pars.dat");
                phmm_aln->phmm = new t_phmm(phmm_pars_fp);
                free(phmm_pars_fp);

                //printf("similarity is %f\n", ml_result->ml_similarity);
                phmm_aln->phmm->set_parameters_by_sim(similarity);
                // phmm->dump_parameters();

                double forward_score = linearAlign->forward(seq1, seq2, phmm_aln->phmm->trans_probs, phmm_aln->phmm->emission_probs, using_prior);
                double backward_score = linearAlign->backward(phmm_aln->phmm->trans_probs, phmm_aln->phmm->emission_probs, using_prior);   
                if (isVerbose) printf("Iter %d Seq-Seq %d-%d Likelihood: %f\n", i_iter, i_seq1, i_seq2, forward_score);

                // threshold
                double threshold = phmm_aln->phmm->get_fam_threshold(similarity);
                
                // alignment probability
                linearAlign->cal_align_prob(forward_score, threshold, beam_aln_results[i_seq1][i_seq2]);

                for (int i = 0; i <= seq1.size(); i++) {
                    for(auto &item : beam_aln_results[i_seq1][i_seq2][i]){
                        int k = item.first;
                        beam_aln_results[i_seq2][i_seq1][k][i].prob = item.second.prob;
                        beam_aln_results[i_seq2][i_seq1][k][i].aln_prob = item.second.aln_prob;
                    }
                }

                delete(phmm_aln->phmm); // This is not needed after the computation any more.
                delete(phmm_aln);
                delete(linearAlign);
            }
        } // i_seq2 loop.
    } // i_seq1 loop.

    // initialize match score
    if(using_prior){
        for(unsigned int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++){
            delete[] match_score_rets[i_seq1];
            match_score_rets[i_seq1] = new match_score[sequences[i_seq1]->numofbases + 1];
        } 
    }

    return(0);
}

double LinearTurboFold::get_match_score(int i_seq1, int i_seq2, int i, int k)
{
    // calculate
    double temp_upstream_info = sqrt(match_score_rets[i_seq1][i-1].upstream * match_score_rets[i_seq2][k-1].upstream);
    double temp_downstream_info = sqrt(match_score_rets[i_seq1][i-1].downstream * match_score_rets[i_seq2][k-1].downstream);
    double temp_unpairing_info = sqrt(match_score_rets[i_seq1][i-1].unpair * match_score_rets[i_seq2][k-1].unpair);

    return (temp_upstream_info + temp_downstream_info) * 1.0 + temp_unpairing_info * 0.8 + 0.5;
}


int LinearTurboFold::GetNumberSequences()
{
    return(sequences.size());
}

void LinearTurboFold::SetProgress(ProgressHandler& Progress) {
    progress = &Progress;
}

LinearTurboFold::~LinearTurboFold()
{
    const int n_seq = sequences.size();

    for(int i_seq = 0; i_seq < n_seq; i_seq++) {
        if (saves[i_seq]!=NULL) delete[] saves[i_seq];
        delete folds[i_seq];
        free(similarities[i_seq]);
    } 
    free(similarities);

    for(unsigned int i_seq1 = 0; i_seq1 < n_seq; i_seq1++)
    {
        for(unsigned int i_seq2 = i_seq1+1; i_seq2 < n_seq; i_seq2++)
        {
            if(i_seq1 != i_seq2)
            {
                t_phmm_aln* phmm_aln = new t_phmm_aln(sequences[i_seq1], sequences[i_seq2]);
                delete(phmm_aln);

                // Free beam_aln_results
                delete[] this->beam_aln_results[i_seq1][i_seq2];
                delete[] this->beam_aln_results[i_seq2][i_seq1];
            }
            else
            {
            }
        } // i_seq2 loop.

        delete[] match_score_rets[i_seq1];
    } // i_seq1 loop.
    delete[] match_score_rets;

    // sequences are used in memory free'ing. Must make sure it is free'ed at the very last step.
    for(int i_seq = 0; i_seq < n_seq; i_seq++)
    {
        if (sequences[i_seq]!=NULL) delete sequences[i_seq];
    }

    if (multiple_sequences!=NULL) delete multiple_sequences;
    if (multiple_alignment!=NULL) delete multiple_alignment;
}


int LinearTurboFold::WriteCt(const int i_seq, const char fp[])
{
    int ret = folds[i_seq]->WriteCt(fp); 

    return 1;
}


// -------------------------------------------------------------

int main(int argc, char** argv){
    // parse config file
    ConfigParser config;
    config.ParseConfig(argc, argv);

    // count runtime if verbose
    struct timeval parse_starttime, parse_endtime; // time
    if (config.verbose) gettimeofday(&parse_starttime, NULL);
    
    // initializing linearturbofold
    LinearTurboFold *LTF;
	LTF = new LinearTurboFold(config.fasta_sequences, config.ckybeam, config.hmmbeam, config.verbose); // TODO: only support fasta input TODO: align beam

    // create the progress monitor.
    TProgressDialog* progress = new TProgressDialog();

    // run iteractions
    LTF->SetProgress( *progress );
    LTF->run_iterations(config); 

    if (config.verbose){
        gettimeofday(&parse_endtime, NULL); 
	    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
	    cout << "LinearTurboFold Time: " << parse_elapsed_time << " seconds." << endl;
    }
    

    return 0;
}
