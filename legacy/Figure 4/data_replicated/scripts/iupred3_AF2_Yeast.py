#!/usr/bin/python3

import sys
import csv
import os
import iupred3_lib
import argparse
from Bio import SeqIO
from os import path


# This function takes in the file path of the tab-separated text file as an argument,
# opens the file, reads it using the csv reader with a delimiter of tab character,
# iterates through each row and appends the string value in the 1st column to the pairs_list.
# Remove 'None'
# Finally, it returns the pairs_list.
def get_list_fasta_AF2(file_path):
    pairs_list_pos = []
    pairs_list_neg = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader) # Skip the first row
        for row in reader:
            pairs_list_pos.append(row[0])
            pairs_list_neg.append(row[2])
                                
        pairs_list_pos = [x for x in pairs_list_pos if x != 'None']
        pairs_list_neg = [x for x in pairs_list_neg if x != 'None']
        pairs_list_pos = [x + '.fasta' for x in pairs_list_pos]
        pairs_list_neg = [x + '.fasta' for x in pairs_list_neg]
    
    return pairs_list_pos, pairs_list_neg


# The function takes in the path to a directory containing FASTA files as an input,
# and reads in all the files in that directory. If the id of the sequence already exists,
# it will skip it.
# The function will return a dictionary where the keys are the ids of the sequences and the
# values are the sequences themselves in the form of SeqRecord objects.
def read_fasta(directory, list_fasta):
    fasta_dict = {}
    for file in list_fasta:
        for record in SeqIO.parse(os.path.join(directory, file), "fasta"):
            if record.id not in fasta_dict:
                fasta_dict[record.id] = record
            else:
                continue
                
    return fasta_dict



PATH = os.path.dirname(os.path.realpath(__file__))
AF2_output_path = '../INTS_outputs/Yeast_avg_ints_ColabFold_AF2Mv2_unpaired+paired.txt'
pair_pos, pair_neg = get_list_fasta_AF2(AF2_output_path)


directory_pos = '../../data/FASTA_files/Yeast/Yeast_Positive_set/'
directory_neg = '../../data/FASTA_files/Yeast/Yeast_Negative_set/'
fasta_dict_pos = read_fasta(directory_pos, pair_pos)
fasta_dict_neg = read_fasta(directory_neg, pair_neg)


for key, value in fasta_dict_pos.items():
    uniprot_id = key
    sequence = value
    s = sequence.seq
    if not path.exists('../iupred3/results_AF2_Yeast_DIP_Pos/{}_iupred3.txt'.format(uniprot_id)):
        iupred2_result = iupred3_lib.iupred(sequence, 'long')
            #print("""# IUPred3 - improved prediction of protein disorder with a focus on specific user applications 
            # Gábor Erdős, Mátyás Pajkos, Zsuzsanna Dosztányi
            # Nucleic Acids Research 2021, Submitted
            # 
            # IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding
            # Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi
            # Nucleic Acids Research 2018;46(W1):W329-W337.
            #
            # Prediction type: {}
            # Smoothing used: {}
            # Prediction output""".format(args.iupred_type, args.smoothing))

        f = open('../iupred3/results_AF2_Yeast_DIP_Pos/{}_iupred3.txt'.format(uniprot_id), 'w')
        for pos, residue in enumerate(sequence):
            f.write('{}\t{}\t{:.4f}\n'.format(pos + 1, residue, iupred2_result[0][pos]))
        f.close()
    else:
        continue

for key, value in fasta_dict_neg.items():
    uniprot_id = key
    sequence = value
    s = sequence.seq
    if not path.exists('../iupred3/results_AF2_Yeast_DIP_Neg/{}_iupred3.txt'.format(uniprot_id)):
        iupred2_result = iupred3_lib.iupred(sequence, 'long')
            #print("""# IUPred3 - improved prediction of protein disorder with a focus on specific user applications 
            # Gábor Erdős, Mátyás Pajkos, Zsuzsanna Dosztányi
            # Nucleic Acids Research 2021, Submitted
            # 
            # IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding
            # Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi
            # Nucleic Acids Research 2018;46(W1):W329-W337.
            #
            # Prediction type: {}
            # Smoothing used: {}
            # Prediction output""".format(args.iupred_type, args.smoothing))

        f = open('../iupred3/results_AF2_Yeast_DIP_Neg/{}_iupred3.txt'.format(uniprot_id), 'w')
        for pos, residue in enumerate(sequence):
            f.write('{}\t{}\t{:.4f}\n'.format(pos + 1, residue, iupred2_result[0][pos]))
        f.close()
    else:
        continue


