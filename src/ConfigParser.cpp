#include "ConfigParser.h"
#include <cmath>

ConfigParser::ConfigParser() {

	// Initialize the calculation type description.
	// calcType = "LinearTurboFold";

	// Initialize the TurboFold iterations.
	turboIterations = 3;

	// Initialize the default alignment output max column number.
	AlnFormat = "Fasta";
        ColumnNumber = 60;

	// Initialize ThreshKnot mode parameters.
		// Initialize the number of ProbKnot iterations.
		pkIterations = 1;

		// Initialize the minimum helix length for a pseudoknot.
		minHelixLength = 1;

		// Initialize the default probability threshold.
		threshold = 0.3;

	// the input sequences are all contained in a single FASTA file.
	fasta_sequences = NULL;
	
	// LinearTurboFold beam size
		// The default value is 0
		// -1 means that the beam search is not applicated to HMM alignment
		hmmbeam = 100;

		// The default value is 0
		// -1 means that the beam search is not applicated to partition function
		ckybeam = 100;
	
	// Verbose mode
		verbose = false;

	// set whether save base pair probabililty 
		saveBpps = false;

	// set whether save parition function
		savePfs = false;
}

void ConfigParser::ParseConfig(int argc, char** argv){
	string inputFasta = argv[1]; 
    string outputDir = argv[2];
    hmmbeam = atoi(argv[3]);
    ckybeam = atoi(argv[4]);
    turboIterations = atoi(argv[5]);
    saveBpps = atoi(argv[6]) == 1;
    savePfs = atoi(argv[7]) == 1;
    verbose = atoi(argv[8]) == 1;
    minHelixLength = atoi(argv[9]);
    pkIterations = atoi(argv[10]);
    threshold = atof(argv[11]);

    fasta_sequences = t_structure::read_multi_seq(inputFasta.c_str(), false); // false as final parameter prevents modification of CT label.
    const unsigned int sequenceCount = fasta_sequences->size(); // from here on, the size of outputCtFiles, outputPfsFiles, etc must match the size of the sequenceFiles vector.
    sequenceFiles.resize(sequenceCount); // (used to determine sequenceCount)

    outputCtFiles.resize(sequenceCount); 
    outputDotFiles.resize(sequenceCount);
    if (saveBpps) outputBppFiles.resize(sequenceCount);
    if (savePfs) outputPfsFiles.resize(sequenceCount);
    const int num_width = log10(sequenceCount) + 1; // width of formatted sequence number (for auto-generated names).
    string tmpFileName = "";       
    for(unsigned int i=0; i<sequenceCount; i++) {
        string name = (*fasta_sequences)[i]->ctlabel;
        trim(name); // remove surrounding whitespace.
        name.insert(0, sfmt("%0*i_", num_width, i+1)); // prefix with formatted sequence number, e.g.  "007_"
        // replace invalid file-name characters. Truncate the name if it is too long, and append the .ct extension.
        
        tmpFileName = createSafeFilename(name, ".ct", true);
        outputCtFiles[i] = outputDir + "/" + tmpFileName;

        tmpFileName = createSafeFilename(name, ".db", true);
        outputDotFiles[i] = outputDir + "/" + tmpFileName;
        
        if (saveBpps) {
            tmpFileName = createSafeFilename(name, ".bpp", true);
            outputBppFiles[i] = outputDir + "/" + tmpFileName;
        }
        
        if (savePfs) {
            tmpFileName = createSafeFilename(name, ".pfs", true);
            outputPfsFiles[i] = outputDir + "/" + tmpFileName;
        }
    }

    OutAln = outputDir + "/" + "output.aln";
}
