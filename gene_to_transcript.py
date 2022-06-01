#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 12:01:59 2021

@author: harry

Creating an updated GTF file. Converting gene ID to transcript ID.
"""
import re
import sys


gtf_path = "/home/harry/Documents/lncRNA_August_21/ST_v9.1_merged.gtf"
gtf_file = open(gtf_path).readlines()
updated_path = "/home/harry/Documents/lncRNA_August_21/ST_v9.1_merged_transcript_updated.gtf"

def converter(input_string, target_pattern, replacement_pattern):
    old_string = re.search(target_pattern, input_string).group(1)
    old_border = target_pattern.split("(.*?)")
    old_pattern = old_border[0] + old_string + old_border[1]

    updated_string = re.search(replacement_pattern, input_string).group(1)
    new_border = target_pattern.split("(.*?)")
    new_pattern = new_border[0] + updated_string + new_border[1]
    output = input_string.replace(old_pattern, new_pattern, 1)    
    
    return(output)

AFM_key_path = "/home/harry/Documents/lncRNA_August_21/afm2afub.key"

# Generate a clean AFM key 
AFM_key_file = open(AFM_key_path).readlines()
AFM_key = {}
for line in AFM_key_file:
    AFM_id = line.split(",")[0]
    AFUB_id = line.split(",")[1]
    if re.search("AFUB_[0-9]{1,6}", AFUB_id):
        clean_AFUB_id = re.search("AFUB_[0-9]{1,6}", AFUB_id).group(0)
    else:
        AFUB_id = AFUB_id.replace("hki_", "")
        clean_AFUB_id = re.search("AFUB_[0-9]{1,6}", AFUB_id).group(0)
    AFM_key[AFM_id] = clean_AFUB_id

with open(updated_path, "w") as f:
    for line in gtf_file:
        newline = line
        if ("transcript_id" in line) and ("gene_id" in line):
            if "MSTRG" in line:
                newline = converter(line, 'gene_id \"(.*?)\"', 'transcript_id \"(.*?)\"' )
                if "AFM" in newline:
                    print("Exists")
                    # Extract AFM code
                    AFM_id = re.search(r'cmp_ref \"(.*?)\"', newline).group(1)
                    newline = newline.replace(AFM_id, AFM_key.get(AFM_id, AFM_id))
                    
            #elif "AFUB" in line:
            #    newline = converter(line, 'transcript_id \"(.*?)\"', 'gene_id \"(.*?)\"' )

        f.write(newline)



