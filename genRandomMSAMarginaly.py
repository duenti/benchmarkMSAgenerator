####Marginally positions####
#1 Á Amides: NQ
#2 Â Aliphatic: GAVLI
#3 & Basic: HKR
#4 Ĥ Hydroxile: STY
#5 Ŝ Sulfur: CM
#6 Ñ Non-Polar: FGVLAIPMW
#7 Ṕ Polar: YSNTQC
#8 ! VeryHydrophobic: LIFWVM
#10 # Hydrophilic: RKNQPD
#11 + Pos. Carregado: KR
#12 _ Neg. Carregado(Acidic): DE
#13 Ũ 60-90: GAS
#14 Õ 108-117: CDPNT
#15 Ĩ 138-154: EVQH
#16 Ẽ 162-174: MILKR
#17 Ã 189-228(Aromáticos): FYW
#18 $ Similarity: ND
#19 % Similarity: QE

##Cada grupo de sequência pertence a uma única subclasse
import sys
import os
from joblib import Parallel, delayed
import multiprocessing
from random import randint, choice
import numpy as np, numpy.random
import random

if len(sys.argv) != 2:
	print("Generate random multiple sequence alignments (Configure global variables before run)\nUsage: python genRandomMSAMarginally.py outputdir\n");
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
prob_marginally_conservation = 0.25 #Probability of a conserved (or locally conserved) position be related to a marginally property. 1-p will be the probability of the position be related to an specific amino acid
prob_gap = 0.05 #Probability of a conserved position be gap

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

#def randomAAgap():
#	aas = "-ACDEFGHIKLMNPQRSTVWY-"
#	return choice(aas)

def randomMarginallyAA(idx):
	if idx == '1':
		p_dist = np.random.dirichlet(np.ones(2),size=1)[0]
		return numpy.random.choice(["N","Q"], p=p_dist)
	elif idx == '2':
		p_dist = np.random.dirichlet(np.ones(5),size=1)[0]
		return numpy.random.choice(["G","A","V","L","I"], p=p_dist)
	elif idx == '3':
		p_dist = np.random.dirichlet(np.ones(3),size=1)[0]
		return numpy.random.choice(["H","K","R"], p=p_dist)
	elif idx == '4':
		p_dist = np.random.dirichlet(np.ones(3),size=1)[0]
		return numpy.random.choice(["S","T","Y"], p=p_dist)
	elif idx == '5':
		p_dist = np.random.dirichlet(np.ones(2),size=1)[0]
		return numpy.random.choice(["C","M"], p=p_dist)
	elif idx == '6':
		p_dist = np.random.dirichlet(np.ones(9),size=1)[0]
		return numpy.random.choice(["F","G","V","L","A","I","P","M","W"], p=p_dist)
	elif idx == '7':
		p_dist = np.random.dirichlet(np.ones(6),size=1)[0]
		return numpy.random.choice(["Y","S","N","T","Q","C"], p=p_dist)
	elif idx == '8':
		p_dist = np.random.dirichlet(np.ones(6),size=1)[0]
		return numpy.random.choice(["L","I","F","W","V","M"], p=p_dist)
	elif idx == '9':
		p_dist = np.random.dirichlet(np.ones(6),size=1)[0]
		return numpy.random.choice(["R","K","N","Q","P","D"], p=p_dist)
	elif idx == '10':
		p_dist = np.random.dirichlet(np.ones(2),size=1)[0]
		return numpy.random.choice(["K","R"], p=p_dist)
	elif idx == '11':
		p_dist = np.random.dirichlet(np.ones(2),size=1)[0]
		return numpy.random.choice(["D","E"], p=p_dist)
	elif idx == '12':
		p_dist = np.random.dirichlet(np.ones(3),size=1)[0]
		return numpy.random.choice(["G","A","S"], p=p_dist)
	elif idx == '13':
		p_dist = np.random.dirichlet(np.ones(5),size=1)[0]
		return numpy.random.choice(["C","D","P","N","T"], p=p_dist)
	elif idx == '14':
		p_dist = np.random.dirichlet(np.ones(4),size=1)[0]
		return numpy.random.choice(["E","V","Q","H"], p=p_dist)
	elif idx == '15':
		p_dist = np.random.dirichlet(np.ones(5),size=1)[0]
		return numpy.random.choice(["M","I","L","K","R"], p=p_dist)
	elif idx == '16':
		p_dist = np.random.dirichlet(np.ones(3),size=1)[0]
		return numpy.random.choice(["F","Y","W"], p=p_dist)
	elif idx == '17':
		p_dist = np.random.dirichlet(np.ones(2),size=1)[0]
		return numpy.random.choice(["N","D"], p=p_dist)
	elif idx == '18':
		p_dist = np.random.dirichlet(np.ones(2),size=1)[0]
		return numpy.random.choice(["Q","E"], p=p_dist)
	else:
		print("ERROR")
		return '-'

def checkAAsubgroup(idx):
	if idx == '1':
		return 'Á'
	elif idx == '2':
		return 'Â'
	elif idx == '3':
		return '&'
	elif idx == '4':
		return 'Ĥ'
	elif idx == '5':
		return 'Ŝ'
	elif idx == '6':
		return 'Ñ'
	elif idx == '7':
		return 'Ṕ'
	elif idx == '8':
		return '!'
	elif idx == '9':
		return '#'
	elif idx == '10':
		return '+'
	elif idx == '11':
		return '_'
	elif idx == '12':
		return 'Ũ'
	elif idx == '13':
		return 'Õ'
	elif idx == '14':
		return 'Ĩ'
	elif idx == '15':
		return 'Ẽ'
	elif idx == '16':
		return '$'
	elif idx == '17':
		return '%'
	elif idx == '18':
		return 'Ã'
	else:
		return idx

def getComposedAAList(idx):
	if idx == '1':
		return ["N","Q"]
	elif idx == '2':
		return ["G","A","V","L","I"]
	elif idx == '3':
		return ["H","K","R"]
	elif idx == '4':
		return ["S","T","Y"]
	elif idx == '5':
		return ["C","M"]
	elif idx == '6':
		return ["F","G","V","L","A","I","P","M","W"]
	elif idx == '7':
		return ["Y","S","N","T","Q","C"]
	elif idx == '8':
		return ["L","I","F","W","V","M"]
	elif idx == '9':
		return ["R","K","N","Q","P","D"]
	elif idx == '10':
		return ["K","R"]
	elif idx == '11':
		return ["D","E"]
	elif idx == '12':
		return ["G","A","S"]
	elif idx == '13':
		return ["C","D","P","N","T"]
	elif idx == '14':
		return ["E","V","Q","H"]
	elif idx == '15':
		return ["M","I","L","K","R"]
	elif idx == '16':
		return ["F","Y","W"]
	elif idx == '17':
		return ["N","D"]
	elif idx == '18':
		return ["Q","E"]
	else:
		print("ERROR")
		return []

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
			if bernoulliTest(prob_gap):
				highcons_residues[pos] = '-'
			if bernoulliTest(prob_marginally_conservation):
				highcons_residues[pos] = str(randint(1, 18))
			else:
				highcons_residues[pos] = randomAA()
			validPositions.remove(pos)

	#Define the main subclass dependent residues
	positions = {}
	for group in main_groups:
		while True:
			countPos = 0
			for pos in validPositions:
				if 	bernoulliTest(prob_subclass_conserved):
					countPos += 1
					validPositions.remove(pos)
					#Generate residues for each subclass in group
					for sb_idx in group:
						if bernoulliTest(prob_gap):
							main_sub_residues[sb_idx][pos] = '-'
						if bernoulliTest(prob_marginally_conservation):
							while True:
								aaidx = str(randint(1, 18))
								aalist = getComposedAAList(aaidx)
								if pos in positions:
									if not bool(set(aalist) & set(positions[pos])):
										main_sub_residues[sb_idx][pos] = aaidx
										positions[pos] = list(set(positions[pos]+aalist))
										break
								else:
									main_sub_residues[sb_idx][pos] = aaidx
									positions[pos] = aalist
									break
						else:
							while True:
								aa = randomAA()
								if pos in positions:
									if aa not in positions[pos]:
										main_sub_residues[sb_idx][pos] = aa
										positions[pos].append(aa)
										break
								else:
									main_sub_residues[sb_idx][pos] = aa
									positions[pos] = [aa]
									break
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
						if bernoulliTest(prob_gap):
							sec_sub_residues[sb_idx][pos] = '-'
						if bernoulliTest(prob_marginally_conservation):
							while True:
								aaidx = str(randint(1, 18))
								aalist = getComposedAAList(aaidx)
								if pos in positions:
									if not bool(set(aalist) & set(positions[pos])):
										sec_sub_residues[sb_idx][pos] = aaidx
										positions[pos] = list(set(positions[pos]+aalist))
										break
								else:
									sec_sub_residues[sb_idx][pos] = aaidx
									positions[pos] = aalist
									break
						else:
							while True:
								aa = randomAA()
								if pos in positions:
									if aa not in positions[pos]:
										sec_sub_residues[sb_idx][pos] = aa
										positions[pos].append(aa)
										break
								else:
									sec_sub_residues[sb_idx][pos] = aa
									positions[pos] = [aa]
									break
			if countPos > 1:
				break

	#print(sec_sub_residues)

	#Make the guide text
	for dic in main_sub_residues:
		count = 0
		for pos, aa in dic.items():
			aa = checkAAsubgroup(aa)
			#print(aa)
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
			aa = checkAAsubgroup(aa)
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
					if aa.isdigit():
						sequence[pos-1] = randomMarginallyAA(aa)
					else:
						sequence[pos-1] = highcons_residues[pos] #High conserved position
				#else: #It's an exception in a conserved position. Verify if is gap or random residue
				#	if bernoulliTest(gap_prob) == 0:
				#		sequence[pos-1] = randomAA()

			#Main subclass dependent
			for pos,aa in main_sub_residues[subi].items(): #Subclass dependent residue
				if bernoulliTest(min_prob_cons): #Test if above minimum frequency of high conservation
					if aa.isdigit():
						sequence[pos-1] = randomMarginallyAA(aa)
					else:
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
							if aa.isdigit():
								seq[pos-1] = randomMarginallyAA(aa)
							else:
								seq[pos-1] = sec_sub_residues[sb_idx][pos]
							sequences[i] = (seqname,seq)
			S_seqs1 = S_seqs2

	#sequences.remove(("",["-"]))
	sequences = [x for x in sequences if x != ("",["-"])]
	sequences.sort(key=lambda tup: tup[0]) 

	#Generate the amino acid probability distribution for each column and make a test. These residues does not have a priori correlation, so they are just random distributed
	for i in range(0,len(sequences[0][1])):
		aatypes = ["-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
		aalist = []
		if i+1 in positions:
			aalist = positions[i+1]
		freeaas = list(set(aatypes)-set(aalist))

		prob_distribution = np.random.dirichlet(np.ones(len(freeaas))/2,size=1)
		for j in range(0,len(sequences)):
			if sequences[j][0] != "":
				if sequences[j][1][i] == "*":
					#Make a test for amino acid acording to the column probability distribution
					aa = np.random.choice(freeaas,p=prob_distribution[0])
					#aa = numpy.random.choice(["-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-"], p=prob_distribution[0])
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