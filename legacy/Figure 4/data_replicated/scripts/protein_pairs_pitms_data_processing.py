#!/usr/bin/python
# -*- coding: utf-8 -*-

# Author: Ilyes Abdelhamid, 2023
# Description: For each model preset/pair mode combination, it calculates 
# average piTM score (PITMS) values for protein pairs in positive and negative sets,
# and saves the results to text files.


from pathlib import Path
from itertools import zip_longest
import numpy as np
import glob
import csv
import re
import sys
import json
import pickle
import itertools

def merge_sets_for_all_variants(path_to_pitms_outputs):

	model_preset = ['AlphaFold2-multimer-v2']
	pair_mode = ['unpaired+paired']

	for m in model_preset:
		for p in pair_mode:
			AFM_outputpath_pos = '..\\..\\data\\Yeast_AFM_output\\Positive_set\\{}\\{}\\'.format(m, p)
			AFM_outputpath_neg = '..\\..\\data\\Yeast_AFM_output\\Negative_set\\{}\\{}\\'.format(m, p)
			listDir = [*glob.glob(AFM_outputpath_pos + '*'), *glob.glob(AFM_outputpath_neg + '*')]
			protein_pairs_pos_name = []
			protein_pairs_neg_name = []
			avg_pitms_pos = []
			avg_pitms_neg = []
			if m == 'AlphaFold2-multimer-v1':
				merged_filename = 'Yeast_avg_pitms_ColabFold_AF2Mv1_{}.txt'.format(p)
			elif m == 'AlphaFold2-multimer-v2':
				merged_filename = 'Yeast_avg_pitms_ColabFold_AF2Mv2_{}.txt'.format(p)

			for d in listDir:
				path_to_stats_all_file = glob.glob(d + '\\stats_all_*.json')
				if 'Positive_set' in d:
					if len(path_to_stats_all_file) == 0:
						raise Exception('No stats_all file identified in {}'.format(d))
					else:
						protein_pairs_pos_name.append(re.search('\\\\([0-9]+_.*_.*).fasta', d).group(1))
						f = open(path_to_stats_all_file[0])
						metric = json.load(f)
						avg_pitms_pos.append(round(np.mean([value for key, value in metric['pitms'].items()]), 3))
						f.close()
				elif 'Negative_set' in d:
					if len(path_to_stats_all_file) == 0:
						raise Exception('No stats_all file identified in {}'.format(d))
					else:
						protein_pairs_neg_name.append(re.search('\\\\([0-9]+_.*_.*).fasta', d).group(1))
						f = open(path_to_stats_all_file[0])
						metric = json.load(f)
						avg_pitms_neg.append(round(np.mean([value for key, value in metric['pitms'].items()]), 3))
						f.close()

			resDict = {'Protein Pairs Pos' : protein_pairs_pos_name, 'Avg PITMS Pos' : avg_pitms_pos, \
			'Protein Pairs Neg' : protein_pairs_neg_name, 'Avg PITMS Neg' : avg_pitms_neg}

			with open(path_to_pitms_outputs + merged_filename, 'w') as f:
				for row in zip_longest(*([key] + value for key, value in resDict.items()), fillvalue = None):
					print(*row, sep = '\t', file=f)



path_to_pitms_outputs = '..\\PITMS_outputs\\'
merge_sets_for_all_variants(path_to_pitms_outputs)

