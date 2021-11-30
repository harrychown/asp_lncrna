#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Swap GTF lines containing AFM number with the resepctive AFUB number

@author: Harry Chown
@date: 15/09/21
"""
def converter(input_string, target_pattern, replacement_pattern):
    old_string = re.search(target_pattern, input_string).group(1)
    old_border = target_pattern.split("(.*?)")
    old_pattern = old_border[0] + old_string + old_border[1]

    updated_string = re.search(replacement_pattern, input_string).group(1)
    new_border = target_pattern.split("(.*?)")
    new_pattern = new_border[0] + updated_string + new_border[1]
    output = input_string.replace(old_pattern, new_pattern, 1)    
    
    return(output)

import re

# Inputs
afm2afub_filename="/home/harry/Documents/lncRNA_sept_21/afm2afub.key"
gtf_filename="/home/harry/Documents/lncRNA_sept_21/ST_v10.3_merged.gtf"

# Read AFUB-AFM key to generate dictionary
afm2afub={}
afm2afub_file=open(afm2afub_filename).readlines()
for line in afm2afub_file:
    line=line.rstrip("\n")
    afm=line.split(",")[0]
    afub=line.split(",")[1]
    afm2afub[afm]=afub
    
# Generate a gene ID-transcript ID dictionary for PCG
gtf_file=open(gtf_filename).readlines()
afub2transcript={}
for line in gtf_file:
    if ('\ttranscript\t' in line) and ('EDP' in line):
        transcript=re.search(r'transcript_id \"(.*?)\"', line).group(1)
        afub=re.search(r'gene_id \"(.*?)\"', line).group(1)
        afub2transcript[afub]=transcript

        


# Replace lines in the GTF
new_gtf_filename="/home/harry/Documents/lncRNA_sept_21/ST_v11_merged.gtf"
new_gtf_file=open(new_gtf_filename, "w")
for line in gtf_file:
    newline=line
    # Replace AFM key with AFUB
    if "AFM" in newline:
        afm=re.search(r'cmp_ref \"(.*?)\"', newline).group(1)
        afub=afm2afub.get(afm, afm+"_transposon")
        newline=newline.replace(afm, afub)
    # Replace MSTRG gene ID with transcript ID
    if ("transcript_id" in newline) and ("gene_id" in newline) and ("MSTRG" in newline):
            newline = converter(newline, 'gene_id \"(.*?)\"', 'transcript_id \"(.*?)\"' )
            old_antisense=re.search(r'cmp_ref \"(.*?)\"', newline)
            if old_antisense:
                old_antisense=old_antisense.group(1)
                if "transcript" not in old_antisense:
                    transcript=afub2transcript.get(old_antisense, old_antisense+"_transposon")
                    newline=newline.replace(afub, transcript)
    new_gtf_file.write(newline)
        
