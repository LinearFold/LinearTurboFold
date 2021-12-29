# LinearTurboFold

This repository contains the C++ source code for the LinearTurboFold project, an end-to-end linear-time algorithm for structural alignment and conserved structure prediction of RNA homologs, which is the first joint-fold-and-align algorithm to scale to full-length SARS-CoV-2 genomes without imposing any constraints on base-pairing distance.

[LinearTurboFold: Linear-time global prediction of conserved structures for RNA homologs with applications to SARS-CoV-2](https://www.pnas.org/content/118/52/e2116269118)

_Proceedings of the National Academy of Sciences_, November 2021.

Sizhen Li, He Zhang, Liang Zhang, Kaibo Liu, Boxiang Liu, David Mathews*, Liang Huang*

\* corresponding authors

# Dependency
gcc 4.8.5 or above; <br>
python2.7 

# Compile
```
Make
```

# Run
LinearTurboFold can be run with:
```
./linearturbofold -i input.fasta -o output_dir [OPTIONS]
```
The input file should be in the FASTA format. Please see [input.fasta](input.fasta) as an example. <br>
Output a multiple sequence alignment and predicted secondary structures in the output directory. 

### OPTIONS
`--it`
The number of iterations (default 3). <br>
`--b1`
The beam size for LinearAlignment (default 100, set 0 for infinite beam). <br>
`--b2`
The beam size for LinearPartition (default 100, set 0 for infinite beam). <br>
`--pf`
Save partition functions for all the sequencs after the last iteration (default False). <br>
`--bpp`
Save base pair probabilities for all the sequencs after the last iteration (default False). <br>
`-v`
Print out alignment, folding and runtime information (default False). <br>
`--th`
Set ThreshKnot threshknot (default 0.3). <br>
`--tkit`
Set ThreshKnot iterations (default 1). <br>
`--tkhl`
Set ThreshKnot minimum helix length (default 3). <br>

### Example
```
./linearturbofold -i input.fasta -o results/ --pf --bpp
100% [==================================================]
3 iterations Done!
Outputing partition functions to files ...
Outputing base pair probabilities to files ...
Outputing multiple sequence alignment to results/output.aln...
Outputing structures to files ...
```

# Evalutation Dataset
We used the [RNAStralign](https://rna.urmc.rochester.edu/publications.html) dataset with known alignments and structures to evaluate LinearTurboFold and benchmarks. 

# SARS-CoV-2 Dataset and Results
The 25 SARS-CoV-2 and SARS-related genomes analyzed in the paper are listed in [samples25.fasta](data/sars-cov-2_data/samples25.fasta). <br>
For further study by experts, 
we provide the whole multiple sequence alignment and predicted structures for all genomes from LinearTurboFold in [sars-cov-2_and_sars-related_25_genomes_msa_structures.txt](sars-cov-2_results/sars-cov-2_and_sars-related_25_genomes_msa_structures.txt). <br>
Each genome corresponds to three lines: sequence name, aligned sequence and aligned structure, respectively. 

