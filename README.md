# benchmarkMSAgenerator

This repository contains two python scripts with the aim of generate random multiple sequence alignments that are going to be used as benchmark for locally conservation analysis methods. One of them considers only pairwise correlations between specific amino acids, the other also includes stereochemical and structural properties.

Before use it, open the script and set the global variables.

## Global Variables

**n_alignments:** Number of alignments that are going to be generated
**min_n_seqs:** Minimum number of sequences in the MSA
**max_n_seqs:** Maximum number of sequences in the MSA
**min_n_positions:** Minimum number of columns
**max_n_positions:** Maximum number of columns
**prob_high_conserved:** Probability of a position be high conserved
**min_prob_cons:** Probability of a conserved position maintain it's conserved amino acid, otherwise an outlier residue are going to be used
**min_n_main_subclass:** Minimum number os type 1 subclasses 
**max_n_main_subclass:** Maximum number of type 1 subclasses
**min_n_sec_subclass:** Minimum number of type 2 subclasses
**max_n_sec_subclass:** Maximum number of type 2 subclasses
**prob_diverge:** Probability of a subclass diverge (Do not share the same positions)
**prob_subclass_conserved:** Probability of a position be highly conserved inside a subclass
**min_seqs_subclass:** Minimum fraction of sequences in a subclass
**prob_marginally_conservation:** Probability of a conserved (or locally conserved) position be related to a marginally property. 1-p will be the probability of the position be related to an specific amino acid
**prob_gap:** = 0.05 #Probability of a conserved position be gap
