import os
import sys
import time
from collections import Counter
import math

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import Consensus_Score, ProfileMatrix, MostProbableKmer

def GreedyMotifSearch(DNA, k, PseudoCounts = False, IgnoreWildCard = True):
    # This function computes the best scoring k-mer motifs from a list of DNA strings, using Greedy Search algorithm
    # This algorithm feels very heuristic and not at all guarantees the best solution
    # PseudoCounts: Boolean variable, determining whether to apply Laplace's Rule of Succession, substituting zeros with small numbers called pseudocounts
    # IgnoreWildCard: Boolean variable, if true, all k-mers with wildcard nucleotides will be ignored when computing best motifs. if False, wildcard is assumed to be any nucleotides
    # input: DNA: a list of DNA strings
    #        k: integer
    # output: bestMotifs: a list of k-mer strings

    # starting with the motifs using the 
    bestMotifs = [dna[0:k] for dna in DNA]
    _ ,bestScore = Consensus_Score(bestMotifs)

    t = len(DNA)
    n = len(DNA[0])
    for i in range(n-k+1):
        Motif = DNA[0][i:i+k]
        Motifs = [Motif, ]
        for j in range(1,t):
            profile = ProfileMatrix(Motifs, PseudoCounts = PseudoCounts)
            highest_kmer , _ = MostProbableKmer(DNA[j], profile, k, IgnoreWildCard = IgnoreWildCard)
            Motifs.append(highest_kmer)
        _ , Score = Consensus_Score(Motifs)
        if Score < bestScore:
            bestMotifs = Motifs
            bestScore = Score
            
    return bestMotifs, bestScore