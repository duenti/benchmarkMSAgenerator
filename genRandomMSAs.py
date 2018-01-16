##Cada grupo de sequência pertence a uma única subclasse
import sys
import os
from joblib import Parallel, delayed
import multiprocessing
from random import randint, choice
import numpy as np, numpy.random
import random

if len(sys.argv) != 2:
	print("Generate random multiple sequence alignments (Configure global variables before run)\nUsage: python genRandomMSAs.py outputdir\n");
	sys.exit()

outputdir= sys.argv[1]

##################Requirements#######################
#Joblib: "pip install joblib" or "https://pythonhosted.org/joblib/installing.html"
##################Global Variables###################
n_alignments = 100 #Number of MSA to generate
min_n_seqs = 1000 #Minimum number of sequences in the MSA
max_n_seqs = 2500 #Maximum number of sequences in the MSA
min_n_positions = 50 #Minimum number of columns
max_n_positions = 100 #Max number of columns
prob_high_conserved = 0.1 #Probability of a position be high conserved
min_prob_cons = 0.9 #Probability of a conserved position maintain it's conserved amino acid, otherwise an exception gap or residue
min_n_main_subclass = 2 
max_n_main_subclass = 5
min_n_sec_subclass = 2
max_n_sec_subclass = 5
prob_diverge = 0.35 #Probability of a subclass diverge (Do not share the same positions)
prob_subclass_conserved = 0.1 #Probability of a position be high conserved inside a subclass
min_seqs_subclass = 0.05 #Minimum fraction of sequences in a subclass (A subclass can be very small and all the sequences can have
#								a high identity, just as occurs in real MSAs)
#####################################################

############Alguns conceitos(REMOVER DEPOIS DO ARTIGO)#############
#-Os alinhamentos gerados não tem a intencao de gerarem sequencias biologicamente corretas, apenas para servir de benchmark para
#metodos de correlacao. Portanto a probabilidade de gerar um triptofano ou uma cisteina, por exemplo, e a mesma.
#-Os alinhamentos gerado nao possuem fragmentos de sequencias, nem alta frequencia de gaps, ja que a maioria dos metodos
#os descartam mesmo

def randomAA():
	aas = "ACDEFGHIKLMNPQRSTVWY"
	return choice(aas)

def randomAAgap():
	aas = "-ACDEFGHIKLMNPQRSTVWY-"
	return choice(aas)

def bernoulliTest(p):
	return np.random.binomial(1, p, 1)

def genMSA(na):
	N_seqs = randint(min_n_seqs,max_n_seqs)
	#print("N = " + str(N_seqs))
	N_cols = randint(min_n_positions,max_n_positions)
	N_main_subs = randint(min_n_main_subclass,max_n_main_subclass)
	N_sec_subs = randint(min_n_sec_subclass,max_n_sec_subclass)
	sequences = [("",["-"])]*N_seqs #Create MSA
	validPositions = list(range(1,N_cols+1))
	highcons_residues = {}
	main_subclasses_probs = []
	sec_subclasses_prob = []
	count = 1
	output = str(na) + "\t"
	main_groups = []
	sec_groups = []
	main_sub_residues = []
	sec_sub_residues = []

	#Verify if min_seqs_subclass and max_n_subclass are compatible
	if min_seqs_subclass * max_n_main_subclass > 0.6:
		print("Your number of subclasses and minimum number of sequences are not compatible. Please reduce the value of max_n_subclass or min_seqs_subclass.")
		return
	if min_seqs_subclass * max_n_sec_subclass > 0.6:
		print("Your number of subclasses and minimum number of sequences are not compatible. Please reduce the value of max_n_subclass or min_seqs_subclass.")
		return

	#Generate a valid sizes distribution for main subclasses
	while(True):
		#print("Trying to generate valid subclasses sizes... " + str(count))
		main_subclasses_probs = np.random.dirichlet(np.ones(N_main_subs),size=1) #Dirichlet distribution
		if np.all(main_subclasses_probs > min_seqs_subclass):
			break
		else:
			count += 1

	#Define main subclasses groups
	count = 0
	temp_groups = []
	for i in range(0,N_main_subs):
		temp_groups.append(i)
		main_sub_residues.append({})
		if 	not bernoulliTest(prob_diverge):
			main_groups.append(temp_groups)
			temp_groups = []
	if len(temp_groups) > 0:
		main_groups.append(temp_groups)
	#print(main_groups)

	#Define secondary subclasses groups
	count = 0
	temp_groups = []
	for i in range(0,N_sec_subs):
		temp_groups.append(i)
		sec_sub_residues.append({})
		if 	not bernoulliTest(prob_diverge):
			sec_groups.append(temp_groups)
			temp_groups = []
	if len(temp_groups) > 0:
		sec_groups.append(temp_groups)
	#print(sec_groups)
	
	#Define high conserved residues
	for pos in validPositions:
		if 	bernoulliTest(prob_high_conserved):
			highcons_residues[pos] = randomAAgap() #It is possible to generate a position conserved by gaps
			validPositions.remove(pos)

	#Define the main subclass dependent residues
	for group in main_groups:
		while True:
			countPos = 0
			for pos in validPositions:
				if 	bernoulliTest(prob_subclass_conserved):
					countPos += 1
					validPositions.remove(pos)
					#Generate residues for each subclass in group
					for sb_idx in group:
						main_sub_residues[sb_idx][pos] = randomAAgap()
			if countPos > 1:
				break

	#Define the secondary subclass dependent residues
	for group in sec_groups:
		while True:
			countPos = 0
			for pos in validPositions:
				if 	bernoulliTest(prob_subclass_conserved):
					countPos += 1
					validPositions.remove(pos)
					#Generate residues for each subclass in group
					for sb_idx in group:
						sec_sub_residues[sb_idx][pos] = randomAAgap()
			if countPos > 1:
				break

	#print(sec_sub_residues)

	#Make the guide text
	for dic in main_sub_residues:
		count = 0
		for pos, aa in dic.items():
			count += 1
			if aa != "-":
				residue = aa + str(pos)
				if output[len(output)-1] == "\t":
					output += residue
				else:
					output += "," + residue
		if(count < 1):
			print("Warning: Small community generated in main subclass")
		output += "\t"

	for dic in sec_sub_residues:
		count = 0
		for pos, aa in dic.items():
			count += 1
			#if aa != "-":
			residue = aa + str(pos)
			if output[len(output)-1] == "\t":
				output += residue
			else:
				output += "," + residue
		if(count < 2):
			print("Warning: Small community generated in secondary subclass")
		output += "\t"


	#Generate the alignment
	last_S_seq = 0
	#print(str(subclasses_probs))
	for subi,p in enumerate(main_subclasses_probs[0]):
		S_seqs = int(N_seqs*p)

		#print(str(last_S_seq) + " " + str(S_seqs + last_S_seq))
		for i in range(last_S_seq,S_seqs+last_S_seq):
			#Create sequence name
			seqname = "sub" + str(subi) + "_" + str(i) + "/1-" + str(N_cols)
			#print(seqname)

			#Create sequence
			sequence = ["*"] * N_cols #Temporary character

			#High conservation positions
			for pos, aa in highcons_residues.items():
				if bernoulliTest(min_prob_cons): #Test if above minimum frequency of high conservation
					sequence[pos-1] = highcons_residues[pos] #High conserved position
				#else: #It's an exception in a conserved position. Verify if is gap or random residue
				#	if bernoulliTest(gap_prob) == 0:
				#		sequence[pos-1] = randomAA()

			#Main subclass dependent
			for pos,aa in main_sub_residues[subi].items(): #Subclass dependent residue
				if bernoulliTest(min_prob_cons): #Test if above minimum frequency of high conservation
					sequence[pos-1] = main_sub_residues[subi][pos]
			#else: #It's an exception in a conserved position. Verify if is gap or random residue
			#	if bernoulliTest(gap_prob) == 0:
			#		sequence[pos-1] = randomAA()
			sequences[i] = (seqname,sequence)
		last_S_seq += S_seqs
	
	#Distribute secondary dependent residues
	#Generate a valid sizes distribution for secondary subclasses
	while(True):
		#print("Trying to generate valid subclasses sizes... " + str(count))
		sec_subclasses_prob = np.random.dirichlet(np.ones(N_sec_subs),size=1) #Dirichlet distribution
		if np.all(sec_subclasses_prob > min_seqs_subclass):
			break
	last_S_seq = 0
	for group in sec_groups:
		random.shuffle(sequences)
		S_seqs1 = 0
		for sb_idx in group:
			prob = sec_subclasses_prob[0][sb_idx]
			S_seqs2 = S_seqs1 + int(N_seqs*prob)
			for i in range(S_seqs1,S_seqs2+1):
				for pos,aa in sec_sub_residues[sb_idx].items(): #Subclass dependent residue
					if bernoulliTest(min_prob_cons): #Test if above minimum frequency of high conservation
						seqname,seq = sequences[i]
						if seq != ["-"]:
							seq[pos-1] = sec_sub_residues[sb_idx][pos]
							sequences[i] = (seqname,seq)
			S_seqs1 = S_seqs2

	#sequences.remove(("",["-"]))
	sequences = [x for x in sequences if x != ("",["-"])]
	sequences.sort(key=lambda tup: tup[0]) 

	#Generate the amino acid probability distribution for each column and make a test. These residues does not have a priori correlation, so they are just random distributed
	for i in range(0,len(sequences[0][1])):
		prob_distribution = np.random.dirichlet(np.ones(22)/2,size=1)
		for j in range(0,len(sequences)):
			if sequences[j][0] != "":
				if sequences[j][1][i] == "*":
					#Make a test for amino acid acording to the column probability distribution
					aa = numpy.random.choice(["-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-"], p=prob_distribution[0])
					sequences[j][1][i] = aa

	# #print(sequences[0])					
	# for name,seq in sequences:
	# 	print(name + "\t" + str("".join(seq)))
	#print(output)

	#Write MSA
	print("Writing alignment " + str(na))
	fw = open(outputdir + str(na) + ".txt","w")

	for line in sequences:
		seqname = line[0]
		seq = "".join(line[1])
		if seqname != "":
			fw.write(seqname + " " + seq + "\n")
		#print(seqname + " " + seq)

	fw.close()

	return output

#num_cores = multiprocessing.cpu_count()
num_cores = 3
results = Parallel(n_jobs=num_cores)(delayed(genMSA)(i) for i in range(n_alignments))

fw = open(outputdir + "guide.txt",'w')

for result in results:
	fw.write(result + "\n")

fw.close()