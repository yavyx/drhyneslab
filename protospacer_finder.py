# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 14:21:34 2017

@author: Javier
"""

def fasta_to_string(fasta_name):
    f = open(fasta_name, "r")
    f_list = f.readlines()
    f.close()
    del f_list[0]
    f_string = ''.join(f_list)
    f_string = f_string.replace("\n", "")
    f_string = f_string.upper()
    return f_string

sequence = fasta_to_string("AF3_contig_8.fasta")
previous_30_bases =             sequence[7527-30:7527]
#Adenine Methylase:             AGTCGAAATCATCGACAACGAAGGCAAAGA
#Glycosylase:                   AGAACTTCGTCTTTTTGAATGATGAGCACT
#Cytosine Deaminase:            TCTTGCGTAAAACTTTAAAGGAATAATGAA
#PPF1 Integrase:                GAATGCGGACCGATTGAAGGTGCGGACCGC 
#PPF1 TranscriptionalRegulator1 GCTGGCCGCAGCGCTCTCCTATTTCGAGAG
#PPF1 TranscriptionalRegulator2 TTCGGGGTGCGCGGGATTGAGAGGTGCCAC
#PPF1 Cytosine Methylase        AGAGGCGGCGCGCGCCGCTGAGGGGGCTCA
#PPF1 Adenine Methyltransferase CCGCAGAATGAAGGTGGGACCGATCATGGC
#AF3 dCMP deaminase             [7527-30:7527]

sequence = sequence[7527:8028]
protospacer_text_name = "AF3_dCMP deaminase_protospacers.fasta"
proto_count = 0
protospacer = ''
temp_sequence = ''

f = open(protospacer_text_name, "w")

for i in range(len(sequence) - 2):
    codon = sequence[i] + sequence[i+1] + sequence[i+2]
    if codon[1] == "G" and codon[2] == "G":
        if i <= 30:
            temp_sequence = previous_30_bases + sequence
            proto_start = i
            proto_count += 1
            protospacer = temp_sequence[proto_start:proto_start+33]
            proto_info = ">Protospacer {} | AF3 | dCMP deaminase:\n{}\n\n".format(proto_count, protospacer)
            f.write(proto_info)
            print(proto_info)
            temp_sequence = ''
        else:
            proto_start = i
            proto_count += 1
            protospacer = sequence[proto_start-30:proto_start+3]
            proto_info = ">Protospacer {} | AF3 | dCMP deaminase:\n{}\n\n".format(proto_count, protospacer)
            f.write(proto_info)
            print(proto_info)

f.close()
