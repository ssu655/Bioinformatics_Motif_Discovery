import os
import sys
import time

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import Neighbors, HammingDistance

def MotifEnumeration(DNA, k, d):
    # This function is an exhaustive search for the k-mer motif, with d mismatches, in a list of DNA strings
    # it is exhaustive as it checks every d-neighbors of every k-mer in the first DNA string, against every other DNA string
    # input: DNA: a list of strings; k, d: integer
    # output: a list of (k,d)-motif
    # This algorithm requires the existance of a (k,d)-motif in each of the DNA strings, it won't work if not.
    Patterns = []
    first_DNA = DNA[0]
    n = len(first_DNA)
    for i in range(n-k+1):
        Pattern = first_DNA[i:i+k]
        Neighborhood = Neighbors(Pattern,d)
        for neighbor in Neighborhood: # loop through all d-neighbors of Pattern
            if neighbor in Patterns: # if the neighbor is already marked as a (k,d)-motif, no need to check if it is one again
                break
            else: # check if neighbor is a (k,d)-motif in the DNA strings
                match = 0
                for j in range(1,len(DNA)): # loop through second to last DNA strings
                    j_DNA = DNA[j]
                    for l in range(n-k+1): # loop through all k-mers in that DNA string
                        if HammingDistance(j_DNA[l:l+k], neighbor)<=d:
                            match = 1
                            break
                        else:
                            match =0
                    if match == 0: # that d-neighbor is not a d-neighbor in j-th DNA string
                        break
                if match == 1:
                    Patterns.append(neighbor)
    m = len(Patterns)
    Patterns = list(set(Patterns))
    mm = len(Patterns)
    return Patterns
