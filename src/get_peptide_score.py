#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 17:49:04 2023

@author: magoncal
"""

import pandas as pd
import numpy as np
import os, shutil, sys, getopt
import time
import re
from itertools import combinations, product, permutations

np.set_printoptions(suppress=True)

SOURCE = "withMUTATIONS/covariance_aa_aa/sources/"

path_in = os.path.join(os.path.dirname(os.getcwd()),'refract-data', SOURCE)
marboxes = ["marRAB", "yba0", "rob", "acnA", "acrAB", "fldB", "fldA", "fpr", "hdeA", "mdtG", "poxB", "purA", "ribA", "slp"]
AA_ORDER1 = ['C','F','W','Y']
AA_ORDER1_IDX = {'C':0,'F':1,'W':2,'Y':3}
AA_ORDER2 = ['A','R','N','D','Q','E','G','H','I','L','K','M','P','S','T','V']
AA_ORDER2_IDX = {'A':0,'R':1,'N':2,'D':3,'Q':4,'E':5,'G':6,'H':7,'I':8,'L':9,'K':10,'M':11,'P':12,'S':13,'T':14,'V':15}
AA_ORDER3 = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
AA_CODE_DICT = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}


def print_result(pept, marboxes, mat, det):
    print(pept)
    for row in range(0, len(mat)):
        if det:
            vals = [marboxes[row]]+mat[row].tolist()
        else:
            vals = [marboxes[row], mat[row,0]]
        print(*vals, sep="; ")


def get_peptide_score(argv):
    
    peptide, file, det = "","","0"
    arg_help = "{0} -p <peptide with 6 amino acids> -i <file with peptides, one per line> -det <0 or 1 for detailed values>".format(argv[0])

    # Dealing with the arguments
    try:
        opts, args = getopt.getopt(argv[1:], "h:p:f:d:", ["help", "peptide=", "input=", "details="])
    except:
        print(arg_help)
        sys.exit(2)
        
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-p", "--peptide"):
            peptide = arg
        elif opt in ("-i", "--input"):
            file = arg
        elif opt in ("-d", "--details"):
            det = arg
    #peptide = 'CARNDA'
    #peptide = 'FVPSRS'
    if peptide!="":
        assert(len(peptide) == 6), "Peptide must have 6 amino acids."
        assert(peptide[0] in AA_ORDER1), "Only four residues are allowed in position 1 (C, F, W, Y)."
        assert(len([item for item in peptide[1:] if item in AA_ORDER2])>0), "Only following residues are NOT allowed in positions 2 to 6: C, F, W, Y."
        peptides = [peptide]
    if file!="":
        assert(os.path.exists(file)), "Provide an existing file with multiple peptides, one per line."
        peptides = np.loadtxt(file, dtype='U6').tolist()
        assert(len(peptides[0])==6), "Provide an existing file with multiple peptides, one per line."
    if det!="":
        assert(det.isnumeric()), "Parameter details must be 0 or 1."
        assert((int(det)== 0)|(int(det)== 1)), "Parameter details must be 0 or 1."
        det = int(det)
               
    marb_c = 0
    mat_scores = np.empty([len(peptides), len(marboxes), 16])
    for marb in marboxes:
        pept_c = 0
        mat_deltaE_1 = np.load(path_in+marb+'_ratios_pos1.npy', allow_pickle = True)
        mat_deltaE_2 = np.load(path_in+marb+'_ratios_posN.npy', allow_pickle = True)
        
        ##################
        # I decided to loop several time in the marbox, to reduce disk access.
        ##################
        for pept in peptides:
            pept = pept.strip()
            # Getting the index for position 1 and the others on the matrix
            peptide1 = AA_ORDER1_IDX[pept[0]]
            peptide_idx = [peptide1]+[AA_ORDER2_IDX[i] for i in list(pept[1:])]
    
            # As we are doing a combination, position res1 has 5 matrices, res2 has 4 matrices and so on
            # We now need the combinations (no repetition) from the 6 positions to define dimension 0 of the 3D matrices
            pairs_pept = list(combinations(peptide_idx, 2))
            # Extracting indexes for matrix 1 (4x16)
            idx1 = np.array([[i]+list(j) for i,j in zip(range(5), pairs_pept[0:5])])
            # Extracting indexes for matrix 2 (16x16)
            idx2 = np.array([[i]+list(j) for i,j in zip(range(10), pairs_pept[5:])])
        
            all_deltaE = np.concatenate((mat_deltaE_1[idx1[:, 0],idx1[:, 1], idx1[:,2]], mat_deltaE_2[idx2[:, 0],idx2[:, 1], idx2[:,2]]))
            sum_deltaE = np.sum(mat_deltaE_1[idx1[:, 0],idx1[:, 1], idx1[:,2]])+np.sum(mat_deltaE_2[idx2[:, 0],idx2[:, 1], idx2[:,2]])
            
            mat_scores[pept_c, marb_c, :] = np.insert(all_deltaE, 0, sum_deltaE)
            pept_c+=1
        marb_c+=1
        
    # Output
    for pept_c in range(0, len(peptides)):
        print_result(peptides[pept_c], marboxes, mat_scores[pept_c], det)
            

if __name__ == "__main__":
    get_peptide_score(sys.argv)