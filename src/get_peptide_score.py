#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 17:49:04 2023

@author: magoncal
"""

import pandas as pd
import numpy as np
import os, shutil
import time
import re
from itertools import combinations, product, permutations

np.set_printoptions(suppress=True)

SOURCE = "withMUTATIONS/covariance_aa_aa/"

path_in = '/home/magoncal/Documents/data/projects/idr_cook/refract_chile/refract-data/'+SOURCE
marboxes = ["marRAB", "yba0", "rob", "acnA", "acrAB", "fldB", "fldA", "fpr", "hdeA", "mdtG", "poxB", "purA", "ribA", "slp"]
AA_ORDER1 = ['C','F','W','Y']
AA_ORDER1_IDX = {'C':0,'F':1,'W':2,'Y':3}
AA_ORDER2 = ['A','R','N','D','Q','E','G','H','I','L','K','M','P','S','T','V']
AA_ORDER2_IDX = {'A':0,'R':1,'N':2,'D':3,'Q':4,'E':5,'G':6,'H':7,'I':8,'L':9,'K':10,'M':11,'P':12,'S':13,'T':14,'V':15}
AA_ORDER3 = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
AA_CODE_DICT = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

#peptide = 'CARNDA'
peptide = 'FVPSRS'
assert(peptide[0] in AA_ORDER1), "Only four residues are allowed in position 1 (C, F, W, Y)."
assert(len([item for item in peptide[1:] if item in AA_ORDER2])>0), "Only following residues are NOT allowed in positions 2 to 6: C, F, W, Y."

# Getting the index for position 1 and the others on the matrix
peptide1 = AA_ORDER1_IDX[peptide[0]]
peptide_idx = [peptide1]+[AA_ORDER2_IDX[i] for i in list(peptide[1:])]

# As we are doing a combination, position res1 has 5 matrices, res2 has 4 matrices and so on
# We now need the combinations (no repetition) from the 6 positions to
# define dimension 0 of the 3D matrices
pairs_pept = list(combinations(peptide_idx, 2))
# Extracting indexes for matrix 1 (4x16)
idx1 = np.array([[i]+list(j) for i,j in zip(range(5), pairs_pept[0:5])])
# Extracting indexes for matrix 2 (16x16)
idx2 = np.array([[i]+list(j) for i,j in zip(range(10), pairs_pept[5:])])

for marb in marboxes:
    mat_deltaE_1 = np.load(path_in+'sources/'+marb+'_ratios_pos1.npy', allow_pickle = True)
    mat_deltaE_2 = np.load(path_in+'sources/'+marb+'_ratios_posN.npy', allow_pickle = True)

    all_deltaE = np.concatenate((mat_deltaE_1[idx1[:, 0],idx1[:, 1], idx1[:,2]], mat_deltaE_2[idx2[:, 0],idx2[:, 1], idx2[:,2]]))
    sum_deltaE = np.sum(mat_deltaE_1[idx1[:, 0],idx1[:, 1], idx1[:,2]])+np.sum(mat_deltaE_2[idx2[:, 0],idx2[:, 1], idx2[:,2]])
    
    print('Marbox: {0}'.format(marb))
    print('Score: {0}\n'.format(str(sum_deltaE)))
    print('All values:\n {0}\n'.format(np.array2string(all_deltaE, separator=', ')))
