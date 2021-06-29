
#include "probknot.h"


// return the pairing probability of the i=j pair, where i<j.
double LinearGetBasePairProb(int i, int j, unordered_map<int, State> *pfscore, unordered_map<int, ExtValue>* ext_info, double viterbi_score){
	if (pfscore[j-1].find(i-1) == pfscore[j-1].end()) return 0.0; // N.B.
    if (ext_info[j][i].score < TO_XLOG(EPSILON)) return 0.0; // -690.776

    // double prob = DIV(PROD(v->f(i,j), v->f(j,i+ct->GetSequenceLength())), PROD(w5[ct->GetSequenceLength()] POWSCALING2(2), ct->extrinsic_info[j][i]));
    double prob = DIV(PROD(pfscore[j-1][i-1].alpha, pfscore[j-1][i-1].beta), PROD(viterbi_score POWSCALING2(2), ext_info[j][i].score));
    if (prob > float(-9.91152)) {
        // if (prob > 0.1) cout << i << " " << j << " " << v->f(i,j) << " " <<  v->f(j,i+ct->GetSequenceLength()) << " " << w5[ct->GetSequenceLength()]<< " " << ct->extrinsic_info[j][i] << " " << prob << endl;
        prob = TO_LINEAR(prob);
        if (prob > 1.0) prob = 1.0; // lisiz: no pairs with prob > 1.0, TBD
        return prob;
    }
    return 0.0;
}

int LinearProbKnotPartition(structure *ct, unordered_map<int, State> *pfscore, unordered_map<int, ExtValue>* ext_info, double viterbi_score, PFPRECISION **probs, PFPRECISION *rowprob ){
	//First determine pair probabilities:
	for (int i=1;i<ct->GetSequenceLength();i++) {
		for (int j=i+minloop+1;j<=ct->GetSequenceLength();j++) {
		
			// probs[j][i] = calculateprobability(i,j,v,w5,ct,data,lfce,mod,scaling,fce);
			probs[j][i] = LinearGetBasePairProb(i, j, pfscore, ext_info, viterbi_score);
			// if (probs[j][i] > 0) cout << "new: "<< i << " " << j << " " << probs[j][i] << endl;

			//also accumulate the best probs for each nucleotide:
			if (probs[j][i]>rowprob[i]) rowprob[i] = probs[j][i];
			if (probs[j][i]>rowprob[j]) rowprob[j] = probs[j][i];

		}

	}
    return 0;
}

int LinearProbKnotAssemble(structure *ct, unordered_map<int, State> *pfscore, unordered_map<int, ExtValue>* ext_info, double viterbi_score, int iterations, int MinHelixLength, double threshold)
{
	PFPRECISION **probs,*rowprob;
	int i,j,iter;
	
	//Add one structure:
	ct->AddStructure();

    //Build a 2-d array for storing pair probabilities, probs, note that the higher index is addressed first...
	probs = new PFPRECISION *[ct->GetSequenceLength()+1];

	//also allocate space for rowprob[i], the highest probability for pairing of nucleotide i
	rowprob = new PFPRECISION [ct->GetSequenceLength()+1];

	for (i=1;i<=ct->GetSequenceLength();i++) {
		probs[i] = new PFPRECISION [i+1];
		
		//Initialize rowprob to zero
		rowprob[i] = 0.0;
	}
    
    //Read the partition function and populate "probs" and "rowprob" arrays with probabilities
    // ProbKnotPartition( v, w5, ct, data, lfce, mod, scaling, fce, probs, rowprob );
	LinearProbKnotPartition(ct, pfscore, ext_info, viterbi_score, probs, rowprob );

    //Calculate maximum expected accuracy structure
    ProbKnotCompute( ct, probs, rowprob, iterations, MinHelixLength, threshold );

	//cleanup memory use:
	for (i=1;i<=ct->GetSequenceLength();i++) delete[] probs[i];
	delete[] probs;

	delete[] rowprob;

	return 0;
}

int ProbKnotCompute( structure *ct, PFPRECISION **probs, PFPRECISION *rowprob, int iterations, int MinHelixLength, double threshold ){

	//now assemble the structure:
	for (int i=1;i<ct->GetSequenceLength();i++) {
		for (int j=i+minloop+1;j<=ct->GetSequenceLength();j++) {
			
			//check all possible pairs
			//take a pair if it has the highest prob for any pair involving i or j
			if (rowprob[i]==probs[j][i]&&rowprob[j]==probs[j][i]&&probs[j][i]>=threshold) {	
				ct->SetPair(i,j);
			}
		}
	}

	//Finally, post-process the structures to remove short helices, if specified:
	if (MinHelixLength>1) {		
        RemoveShortHelices(ct, MinHelixLength, 1);
	}
    return 0;
}

//Remove short helices, allowed stacks across single bulges
//Implemented by Stanislav Bellaousov.
void RemoveShortHelices(structure *ct, int MinHelixLength, int StructureNumber) {
	int pairs,i,j;	
	
	//Checking for helixes smaller then argv[5]
	for (i=1;i<=ct->GetSequenceLength();i++) {
	  //	  if(ct->basepr[StructureNumber][i]!=0){
	  if (ct->GetPair(i,StructureNumber)>i){
	    j=ct->GetPair(i,StructureNumber);
	    pairs=1;
	    while (ct->GetPair(i+1,StructureNumber)==j-1||ct->GetPair(i+2,StructureNumber)==j-1||ct->GetPair(i+1,StructureNumber)==j-2){
	      if (ct->GetPair(i+1,StructureNumber)==j-1){
			i++;
			j--;
			pairs++;
	      }
	      else if (ct->GetPair(i+2,StructureNumber)==j-1){
		if (ct->GetPair(i+1,StructureNumber)!=0){
		    ct->RemovePair(ct->GetPair(i+1,StructureNumber),StructureNumber);
		    ct->RemovePair(i+1,StructureNumber);
		  }
			i=i+2;
			j--;
			pairs++;
	      }
	      else {
			i++;
			j=j-2;
			pairs++;
	      }
	    }


	   
	    //Deleting helixes smaller then MinHelixLength

	    if (pairs<MinHelixLength){
			//ct->RemovePair(ct->GetPair(i,StructureNumber),StructureNumber);
			ct->RemovePair(i,StructureNumber);
			
			if(i>=3){
			  while (ct->GetPair(i-1,StructureNumber)==j+1||ct->GetPair(i-2,StructureNumber)==j+1||ct->GetPair(i-1,StructureNumber)==j+2){
			    if (ct->GetPair(i-1,StructureNumber)==j+1){
			      
					ct->RemovePair(ct->GetPair(i-1,StructureNumber),StructureNumber);
					ct->RemovePair(i-1,StructureNumber);
					
					i--;
					j++;
			    }
			    else if (ct->GetPair(i-2,StructureNumber)==j+1){
					ct->RemovePair(ct->GetPair(i-2,StructureNumber),StructureNumber);
					ct->RemovePair(i-2,StructureNumber);
			      
					i=i-2;
					j++;
			    }
			    else {
					ct->RemovePair(ct->GetPair(i-1,StructureNumber),StructureNumber);
					ct->RemovePair(i-1,StructureNumber);
					
			      
					i--;
					j=j+2;
			      
			    }
			    
			  }
			}
			else if(i==2){
			  while (ct->GetPair(i-1,StructureNumber)==j+1||ct->GetPair(i-1,StructureNumber)==j+2){
			    if (ct->GetPair(i-1,StructureNumber)==j+1){
					ct->RemovePair(ct->GetPair(i-1,StructureNumber),StructureNumber);
					ct->RemovePair(i-1,StructureNumber);
			      
					i--;
					j++;
			    }
			    else {
					ct->RemovePair(ct->GetPair(i-1,StructureNumber),StructureNumber);
					ct->RemovePair(i-1,StructureNumber);
			      
					i--;
					j=j+2;
			      
			    }
			    
			  }
			}
	    }
	  }
	}
}
