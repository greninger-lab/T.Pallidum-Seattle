#!/usr/bin/env python3

import subprocess
import argparse
import sys
import os
import regex
from itertools import chain
import numpy as np
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
import pandas as pd

#This script is ONLY to mask for phylogeny. Does not include the interim "lightly masked" for upload to NCBI.
#intra-rRNA tRNA, arp repeats, tp0470 repeats, tprA, tprB, tprC, tprD, tprEFG (inc. duf2715), tprH, tprIJ (inc duf 2715), tprK (whole gene), tprL

# Matches a read to a specified string of nucleotides, "primer".
def fuzzy_match(read_seq, primer,num_mismatches):
	# Finds exact match first.
	exact_match = regex.search(primer,read_seq)
	if exact_match:
		return exact_match[0].rstrip()
	# If can't find exact match, searches for best match with less than specified number of mismatches.
	else:
		if(num_mismatches==1):
			fuzzy_match = regex.search(r"(?b)("+primer + "){s<=1}", read_seq)
		else:
			fuzzy_match = regex.search(r"(?b)("+primer + "){s<=3}", read_seq)
		if fuzzy_match:
			return fuzzy_match.group()
		else:
			return 0
			#if(len(fuzzy_match)>1):
				#print("!!!!!!!!! Multiple matches found!!!!!!!!!")
				#print(fuzzy_match.group())
				#print(fuzzy_match[0])
				#print(fuzzy_match[1])
				#return 0
			#return fuzzy_match[0]

# Masks repetitive/variable regions based on upstream and downstream elements.
def mask_variable_region(fasta_seq,upstream,downstream,gene_name,num_mismatches):
	# Finds the upstream portion of the gene.
	gene_upstream = fuzzy_match(fasta_seq,upstream,num_mismatches)

	#print(upstream)
	#print(fasta_seq)
	#print(gene_upstream)
	#print(num_mismatches)

	if (gene_upstream == 0):
		#print(gene_name,"not found.")
		return fasta_seq

	# Finds the start position of the gene of interest.
	gene_start = str.index(fasta_seq, gene_upstream) + len(gene_upstream)

	# Searches for the downstream portion of the gene starting from the start of the gene to a certain number of bases away, depending on the gene,
	# to further limit multiple matches.
	if(gene_name == "arp gene"):
		fuzzy_gene_end = gene_start + 3000
	elif(gene_name == "tp0470 gene"):
		fuzzy_gene_end = gene_start + 3000
	else:
		fuzzy_gene_end = gene_start + 10000
	gene_downstream = fuzzy_match(fasta_seq[gene_start:fuzzy_gene_end],downstream,num_mismatches)

	if (gene_downstream == 0):
		#print(gene_name,"not found.")
		return fasta_seq

	# Finds the end position of the gene
	gene_end = str.index(fasta_seq,gene_downstream)

	# Grabs the sequence between gene start and end positions
	gene_seq = fasta_seq[gene_start:gene_end]

	#print(gene_name,"found.")#: ",gene_seq,". Masking...")

	# Replace the found genes with the same length of Ns
	return fasta_seq.replace(gene_seq,'N'*len(gene_seq))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Masking variable regions in for Tpallidum WGS.')
	parser.add_argument('-f', '--fasta', help='Fasta file to mask variable regions for.')

	# Checks for argument sanity.
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(1)

	fasta_file = args.fasta

	masked_sequence = ""

	sample_name = fasta_file.split(".fasta")[0]
	#print("Masking genes in ",fasta_file,"...")
	masked_fasta_file = open(sample_name+"_masked_phylo.fasta","w+")

	for line_num, line in enumerate(open(fasta_file)):
		# Writes fasta header
		if(line_num==0):
			masked_fasta_file.write(line)
		# Reads fasta sequence and masks variable/repetitive regions
		else:

			masked_sequence = masked_sequence + line

	masked_sequence = masked_sequence.replace("\n","")

	# intra-rrna trna1 gene
	seq_to_replace = masked_sequence[230000:235000]
	replaced_gene = mask_variable_region(seq_to_replace,"TCTCCCCTTCCCTTTTGAAAA","NNNNNNNNNNNNNNNNNNCTTTG","intra-rrna trna1",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	# intra-rrna trna2 gene
	seq_to_replace = masked_sequence[279000:285000]
	#replaced_gene = mask_variable_region(seq_to_replace,"TCTCCCCTTCCCTTTTGAAAA","CTATTATTCTTTATGTCCCTT","intra-rrna trna2",1)
	replaced_gene = mask_variable_region(seq_to_replace,"TCTCCCCTTCCCTTTTGAAAA","GGGAAAGCCCACTATTATTC","intra-rrna trna2",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)


	#remask
	# arp gene
	#seq_to_replace = masked_sequence[457255:469137]
	seq_to_replace = masked_sequence[462000:467000]
	#replaced_gene = mask_variable_region(seq_to_replace,"AGCTGGGAGGAGTCTT","TCAGGTCGGCCCTCCTC","arp gene",1)
	#Make end of gene manual
	replaced_gene = mask_variable_region(seq_to_replace,"AGCTGGGAGGAGTCTT","TCAGGTCGGCCCTCCTC","arp_gene",1)
	#replaced_gene = mask_variable_region(seq_to_replace,"TTTGGTTTCCCCTTTGTCTC","AGGTCGCTTCTCAGCATACG","arp gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)



	##tp0470 gene
	seq_to_replace = masked_sequence[493859:504150]
	replaced_gene = mask_variable_region(seq_to_replace,"GCGCGCTTGAGAGCTTCAAA","GCGCTGCAGCCNCTGNTC","tp0470 gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	#tprA
	seq_to_replace = masked_sequence[7000:12000]
	replaced_gene = mask_variable_region(seq_to_replace,"AATGTCGCCGTCTCCCTGGG","NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAATG","tprA gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	#tprB
	seq_to_replace = masked_sequence[9000:13500]
	replaced_gene = mask_variable_region(seq_to_replace,"AGGGTCATAATGGCTGCCCT","GGGAGCCTTCAATAAGGGGG","tprB gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	#tprC
	seq_to_replace = masked_sequence[132000:139000]
	replaced_gene = mask_variable_region(seq_to_replace,"GAATCGAACAGCACGCGTTT","TCACTCCCCCTCCTCTCACT","tprC gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	##tprD
	seq_to_replace = masked_sequence[150000:157000]
	replaced_gene = mask_variable_region(seq_to_replace,"CACGGGGAGTGCACGCGCTT","TCACTCCCCCTCCTCTCACT","tprD gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	##tprEFG locus - starts at end of tp0312
	seq_to_replace = masked_sequence[328000:340000]
	replaced_gene = mask_variable_region(seq_to_replace,"GCAGGTGCAAGGAGTGGTGAGT","NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTCGATAAAGGGG","tprEFG locus",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)


	#remask
	##tprH
	#seq_to_replace = masked_sequence[660000:667000]
	seq_to_replace = masked_sequence[660000:668500]
	replaced_gene = mask_variable_region(seq_to_replace,"CGGACTTTTAGGATGGAGCC","AACGCCAAGGTAGGACATAC","tprH gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	#tprIJ locus
	seq_to_replace = masked_sequence[669000:679000]
	#replaced_gene = mask_variable_region(seq_to_replace,"GTTCCTCCCGGCGCCTCCC","ACCCCATGCTACCTCACCCC","tprIJ locus",1)
	replaced_gene = mask_variable_region(seq_to_replace,"GTTCCTCCCGGCGCCTCCC","TCTACCTCNCCCCCCCCCNCCCNNCT","tprIJ locus",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)


	#remask
	#tprK
	seq_to_replace = masked_sequence[973000:981500]
	replaced_gene = mask_variable_region(seq_to_replace,"AAGAAAAGAACCATACATCC","TCAGAATCCGGAACTGCGTT","tprK gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	#remask
	#tprL
	#seq_to_replace = masked_sequence[1123000:1130000]
	seq_to_replace = masked_sequence[129000:1320000]
	#replaced_gene = mask_variable_region(seq_to_replace,"GCGGCGGGGCGCGCTCAGGC","CGGGGCGCCGCCTGCGGACT","tprL gene",1)
	replaced_gene = mask_variable_region(seq_to_replace,"GGGGGGGGGTGTTGTAAAA","CGGGGCGCCGCCTGCGGACT","tprL gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	#remask
	#Messed up allignment
	#seq_to_replace = masked_sequence[1123000:1130000]
	seq_to_replace = masked_sequence[940000:950000]
	#replaced_gene = mask_variable_region(seq_to_replace,"GCGGCGGGGCGCGCTCAGGC","CGGGGCGCCGCCTGCGGACT","tprL gene",1)
	replaced_gene = mask_variable_region(seq_to_replace,"GGTCACCACCACGTGCTTCT","GAATCGCGGTAGCCAAC","Region gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	masked_fasta_file.write(masked_sequence)
