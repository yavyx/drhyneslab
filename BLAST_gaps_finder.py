# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 16:17:59 2017

@author: Javier
"""

import os
import csv
from Bio.Blast import NCBIXML

def get_gaps(blast_records_list):
    gaps = [[] for a in blast_records_list]
    for i, blast_record in enumerate(blast_records_list):
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                gaps[i].append(hsp.query_start + 1)
                gaps[i].append(hsp.query_end - 1)
    return gaps

def order_start_end(start, end):
    if start > end:
        return end, start
    else:
        return start, end

#Parse BLAST results
os.chdir('C:\\Users\\Javier\\OneDrive\\Dr. Hynes Lab (Calgary)\\BLAST')
result_handle = open("P10VF_vs_AF3_BLAST.xml")
blast_records = NCBIXML.parse(result_handle)
output_text_name = "AF3_target_genes.txt"
f = open(output_text_name, "w")

min_gap_length = 200
overlap_threshold = 100

blast_records_list = []

#Make a list to be able to iterate over many times (blast_records is a one-time use generator)
for blast_record in blast_records:
   blast_records_list.append(blast_record)
   
#Initialize lists and variables
gap_positions = []
gap_lengths = [[] for a in blast_records_list]
long_gaps = [[] for a in blast_records_list]
sorted_genome = [[] for a in blast_records_list]  
target_genes = []
genome_info_list = []
gene_info = ""

#Get gaps from BLAST results
gaps_list = get_gaps(blast_records_list)
    
#Sort gaps and remove 1st and last values
for gap in gaps_list:
    gap.sort()
    if gap != []:
        del gap[0]
        del gap[-1]
    
#Get tuples of every gap start-end positions
for contig in gaps_list:
    gap_positions.append(list(zip(contig[0::2], contig[1::2])))
    
#Get gap lengths for each gap
for i, contig in enumerate(gap_positions):
    for n in contig:
        gap_lengths[i].append(n[1] - n[0]) 

#Get new list with gaps bigger than the threshold
for j, contig in enumerate(gap_lengths):
    for k, length in enumerate(contig):
        if length >= min_gap_length:
            long_gaps[j].append(gap_positions[j][k])

#Parse annotated gene info
with open("AF3_annotated.csv") as genome:
    readCSV = csv.reader(genome, delimiter = ",")
    for gene in readCSV:
        genome_info_list.append(gene)

#Sort annotated genes into lists by contig (to coincide with the long_gaps list)
for gene in genome_info_list:
    for contig in range(len(long_gaps)):
        if gene[0] == "".join(["AF3_c", str(contig)]):
            sorted_genome[contig].append(gene)
            
#Remove first empty list and create another one at the end
del sorted_genome[0]
sorted_genome.append([])

f.write("Parameters: \n")
f.write("Gap length minimum: {} \n".format(min_gap_length))
f.write("Overlap threshold: {} \n\n".format(overlap_threshold))

#Get the genes that overlap with the gaps, according to the overlap threshold
for l, contig in enumerate(sorted_genome):
    for gene in contig:
        start = int(gene[4])
        end = int(gene[5])
        start, end = order_start_end(start, end)
        for m, contig2 in enumerate(long_gaps):
            if l == m:
                for gap in contig2:
                    gap_range = set(range(gap[0], gap[1]))
                    gene_range = set(range(start, end))
                    if (len(gap_range & gene_range) > overlap_threshold) and (gene not in target_genes):
                        gene_info = str(gene[0]) + " " + str(gene[1]) + " " + \
                                    str(gene[3]) + " " + str(gene[4]) + " " + \
                                    str(gene[5]) + " " + str(gene[7])
                        print(gene_info)
                        f.write(gene_info)
                        f.write("\n")
                        target_genes.append(gene)
    
result_handle.close()
f.close()
