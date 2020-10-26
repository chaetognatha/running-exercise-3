# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 09:44:39 2020

@author: ortis
"""

#!/usr/bin/python3
'''
    Title: similarity.py
    Date: 2020-10-23
    Author: Julia Ortis Sunyer
    Description:
        This program will create a similarity matrix using the multiple sequence
        alignments created in the first part of the exercise. The program will create one
        euclidian distance matrix for the identity score and one for the alignment score. 
    List of functions:
        None
    List of non-standard modules:
        None
    Procedure:
        
    Usage:
        score.py input_alignment_fasta_file output_file
'''   

import sys

scoring_file = sys.argv[1]
idn_output = sys.argv[2]
align_output = sys.argv[3]

two_ids_identity_dict = {} #a dictionary for the identity scores with two ids
two_ids_alignment_dict = {} #a dictionary for the alignment scores with two ids
id_set = set() #a set with the ids so that I do not have repeated ids

with open(scoring_file, 'r') as trial: #open part 1 document
    trial.readline() #takes out the first line (the header) found in the file
    for line in trial: #For loop to iterate through the msa file
        line = line.strip().split('\t') #assigns to a varibale line each line it reads stripped from new line character splitted by tabs
        #print(line)
        id_matches = '{}-{}'.format(line[0], line[1]) #creates a variable id_matches and assigns it the two ids
        #print(id_matches)
        two_ids_identity_dict[id_matches] = float(line[2].strip('%')) #adds the identity float number without the % into the right key in two ids identity dictionary
        two_ids_alignment_dict[id_matches] = float(line[3]) #adds the alignment score to the right key in two ids alignment dictionary
        #print(two_ids_identity_dict)
        #print(two_ids_alignment_dict)
        id_set.add(line[0]) #adds the value found in line index 0 to the id set
        id_set.add(line[1]) #adds the value found in line index 1 to the id set
#print(id_set)

multiple_id_list = two_ids_identity_dict.keys() #I create a multiple id list variable to which I assign the keys from two ids identity dictionary (the keys from the two ids alignment dictionary are the same)
#print(multiple_id_list)
multiple_id_list = sorted(multiple_id_list) #I sort the multiple id list
#print(multiple_id_list)
id_list = list(id_set) #I create a list from the id set
id_list = sorted(id_list) #I sort the id list
#print(id_list)
single_id_identity_dict = {} #I create a single id identity dictionary
single_id_alignment_dict = {} #I create a single id alignment dictionary

for ID in id_list: #For loop to iterate through the id list
    single_id_identity_dict[ID] = [] #each id is added to the key in the single id identity dictionary
    single_id_alignment_dict[ID] = [] #each id is added to the key in the single id alignment dictionary
    #print(single_id_dict)
    for key in multiple_id_list: #I iterate through the items in multiple id list
        if ID in key: #if the id is found in the keys found in the multiple id list
            identity_score = two_ids_identity_dict[key] #variable identity score created and assigned the key for two ids identity dictionary
            alignment_score = two_ids_alignment_dict[key] #variable alignment score created and assigned the key for the two ids alignment dictionary
            #print(identity_score)
            #print(alignment_score)
            single_id_identity_dict[ID].append(identity_score) #the identity score is appended to the correct key in the single id identity dictionary
            single_id_alignment_dict[ID].append(alignment_score) #the alignment score is appended to the correct key in the single id alignmetn dictionary
print(single_id_identity_dict)
print(single_id_alignment_dict)

length_list = len(single_id_identity_dict[id_list[0]]) #a lenght_list variable is created which contains the length of the values assigned to the first key in the single id identity dictionary
#print(length_list)
for i in range(length_list + 1): #for i in the range of the length list plus one. I do this to add 0s to the diagonal of the matrix where a sequence is compared to itself
    #print(i)
    single_id_identity_dict[id_list[i]][i:i] = [0] #I do what explained before for the identity dictionary (add 0s in the right places)
    single_id_alignment_dict[id_list[i]][i:i] = [0] #I do what explained before for the alignment dictionary (add 0s in the right places)

for values in single_id_identity_dict.values(): #I check the values in the single id identity dictionary to see if they are ok
    print(values)

for values in single_id_alignment_dict.values(): #I check the values in the single id alignment dictionary to see if they are ok
    print(values)
    
import pandas as pd #I import pandas to make the matrices
DF_var_identity = pd.DataFrame.from_dict(single_id_identity_dict).T #I assign to the variable DF_var_identity the dataframe from the single id identity dictionary using pandas
DF_var_identity.columns = id_list #I assigned to the names of the columns of the matrix the id list
print(DF_var_identity)
DF_var_alignment = pd.DataFrame.from_dict(single_id_alignment_dict).T #I assign to the variable DF_var_alignment the dataframe from the single id alignment dictionary using pandas
DF_var_alignment.columns = id_list #I assigned to the names of the columns of the matrix the id list
print(DF_var_alignment)

from scipy.spatial.distance import pdist, squareform #I use this to calculate the euclidean distances on the values in my matrices
identity_dist = pd.DataFrame(squareform(pdist(DF_var_identity.iloc[:, 1:])), columns=id_list, index=id_list)
print(identity_dist)
alignment_dist = pd.DataFrame(squareform(pdist(DF_var_alignment.iloc[:, 1:])), columns=id_list, index=id_list)
print(alignment_dist)

with open(idn_output, 'w') as identity_output: #I print the identity distance matrix output to the correct file
    print(identity_dist, file=identity_output)

with open(align_output, 'w') as alignment_output: #I print the alignment distance matrix output to the correct file
    print(alignment_dist, file=alignment_output)


