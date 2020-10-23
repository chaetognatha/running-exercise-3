# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 09:44:39 2020

@author: ortis
"""

#!/usr/bin/python3
'''
    Title: score.py
    Date: 2020-10-15
    Author: Julia Ortis Sunyer
    Description:
        This program will first create a similarity matrix using the multiple sequence
        alignments created in the first part of the exercise. The program will create one
        matrix for the identity score and one for the alignment score. Moreover, this code
        will allow the user to input an individual name to then read the table and output
        the most similar individual.
    List of functions:
        
    List of non-standard modules:
        I am using 3 modules:
            1. import sys, to be able to get the files in the terminal.
            2. import pathlib, to check if an output file already exists.
            3. import os, to check if the fna (input) file is empty.
    Procedure:
        
    Usage:
        score.py input_alignment_fasta_file output_file
'''   

two_ids_identity_dict = {}
id_set = set()

with open('try.fna', 'r') as trial:
    trial.readline()
    for line in trial: #For loop to iterate through the msa file
        line = line.strip().split('\t')
        #print(line)
        id_matches = '{}-{}'.format(line[0], line[1])
        #print(id_matches)
        two_ids_identity_dict[id_matches] = float(line[2].strip('%'))
        #print(two_ids_identity_dict)
        id_set.add(line[0])
        id_set.add(line[1])
#print(id_set)

multiple_id_list = two_ids_identity_dict.keys()
#print(multiple_id_list)
multiple_id_list = sorted(multiple_id_list)
#print(multiple_id_list)
id_list = list(id_set)
id_list = sorted(id_list)
#print(id_list)
single_id_dict = {}

for ID in id_list:
    single_id_dict[ID] = []
    #print(single_id_dict)
    for key in multiple_id_list:
        if ID in key:
            identity_score = two_ids_identity_dict[key]
            #print(identity_score)
            single_id_dict[ID].append(identity_score)
#print(single_id_dict)

length_list = len(single_id_dict[id_list[0]])
#print(length_list)
for i in range(length_list + 1):
    #print(i)
    single_id_dict[id_list[i]][i:i] = [0]

for values in single_id_dict.values():
    print(values)

import pandas as pd
DF_var = pd.DataFrame.from_dict(single_id_dict).T
DF_var.columns = id_list
print(DF_var)