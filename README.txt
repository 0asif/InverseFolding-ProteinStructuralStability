Author: Muhammad Asif Ali
Email : maa32@illinois.edu

Supervisor: Dr. Gustavo Caetano-Anolles

Course: CPSC499
Instructor: Dr. Ersoz Elhan

Title: AI Guided Mutation Prediction and Vaccine Development Using Inverse Folding

Comments: This project utilized ProteinMPNN as a starting point to analyze Inverse Folding's capability to predict mutation sites and stable regions 
using the SARS-CoV2-2 Spike Protein as a test case

The structure for the folders/files are:

	-data:
		-input:
			spike_alphafold3/: Directory containing output of AlphaFold3 for S-Protein
			exposed_aa.txt : List of exposed amino acid positions
			helices.txt : Amino acid positions that are part of alpha helices
			IUPredScores.txt : Result output from IUPred3 (https://iupred3.elte.hu/)
			PMPNN_spike_af3.fa : Output from ProteinMPNN
			sheets.txt : Amino acid positions that are part of beta sheets
			spike.fasta : Fasta sequence of SARS-CoV-2 Spike protein (EPI_ISL_402124) termini removed
			spike_alphafold3 : Folder containing the output from AlphaFold3 used for InverseFolding
			spike_full_length.fasta : Fasta sequence of SARS-CoV-2 Spike protein (EPI_ISL_402124) termini included
			spike_mutations.txt: List of mutations on the spike protein during the covid-19 pandemic (https://www.mdpi.com/article/10.3390/biology13030134/s1 (Table 1))

		-output:
			helices_list.txt : List of positions involved in helices (easy parsing)
			mutation_counts.csv : Count table recording the number of times each of the 20 amino acids are recorded at a position
			mutation_frequencies.csv : Frequency table recording the number of times each of the 20 amino acids are recorded at a position
			mutation_frequencies_trasposed.csv : Frequency table transposed for analysis ahead
			sheets_list.txt : List of positions involved in sheets (easy parsing)

	-figures:
		entropy_levels_mutations.tiff
		entropy_per_position.tiff
		softmax_probabilities.tiff

	-scripts:
		-extract_SS.py : Converts the Chimera generated file to a simple text list to record secondary structure
		-get_aa_frequencies.py : This will read the PMPNN output and generate the count and frequency tables
		-get_entropy_probs.py : This will generate the entropy and probability visuals
		-get_heatmap.py : This will generate the main heatmap visualizing the positions retained at various entropy levels
		-run_PMPNN.sh : You will need to run this separately outside of the structure of this project


