# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 12:04:35 2017

@author: Javier
"""

import os
from Bio.Blast import NCBIXML

os.chdir('C:\\Users\\Javier\\OneDrive\\Dr. Hynes Lab (Calgary)\\BLAST')
result_handle = open("AF3_dCMPdeaminase_vs_3841_BLAST.xml")
blast_records = NCBIXML.parse(result_handle)
output_text_name = "AF3_dCMPdeaminase_target_protospacers.txt"
f = open(output_text_name, "w")

allowed_matches = 3
mismatches = 2
PAM_cutoff = 28
proto_count = 0
match_count = 0
matches_per_protospacer = {}
blast_records_list = []

for blast_record in blast_records:
   blast_records_list.append(blast_record)

for blast_record in blast_records_list:
    proto_count += 1
    match_count = 0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            match_count += 1 
            matches_per_protospacer[proto_count] = match_count

proto_count = 0

print("Parameters:")
print("Allowed number of matches: {}".format(allowed_matches))
print("Minimum mismatches: {}".format(mismatches))
print("Proximity to PAM site (bases 31-33): {} \n\n".format(PAM_cutoff))

f.write("Parameters: \n")
f.write("Allowed number of matches: {} \n".format(allowed_matches))
f.write("Minimum mismatches: {} \n".format(mismatches))
f.write("Proximity to PAM site (bases 31-33): {} \n\n".format(PAM_cutoff))

for blast_record in blast_records_list:    
    proto_count += 1
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if matches_per_protospacer[proto_count] <= allowed_matches:
                if (hsp.query_end <= PAM_cutoff) or ((hsp.align_length - hsp.identities) >= mismatches):
                    print("Protospacer {}".format(proto_count))
                    print("Score:", hsp.score)
                    print("Bits:", hsp.bits)
                    print("Expect:", hsp.expect)
                    #print(hsp.num_alignments)
                    print("Identity:", hsp.identities, "/", hsp.align_length)
                    #print(hsp.gaps)
                    #print(hsp.strand)
                    #print(hsp.frame) 
                    #print(hsp.query)
                    print("Start-End:", hsp.query_start, "-", hsp.query_end)
                    #print(hsp.match)
                    print("\n")
                        
                    f.write("Protospacer {} \n".format(proto_count))
                    f.write("Score: {} \n".format(hsp.score))
                    f.write("Bits: {} \n".format(hsp.bits))
                    f.write("Expect: {} \n".format(hsp.expect))
                    f.write("Identity: {}/{} \n".format(hsp.identities, hsp.align_length))
                    f.write("Gaps: {}/{} \n".format(hsp.gaps, hsp.align_length))
                    f.write("Start-End: {}-{} \n".format(hsp.query_start, hsp.query_end))
                    f.write("\n")               
                else: 
                    print("Protospacer {}".format(proto_count))
                    f.write("Protospacer {} \n".format(proto_count))
                    print("MATCH CLOSE TO PAM SITE \n\n")
                    f.write("MATCH CLOSE TO PAM SITE \n\n")
                    break
            
f.close()
result_handle.close()
            
