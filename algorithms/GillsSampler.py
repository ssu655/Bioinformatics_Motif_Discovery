import os
import sys
import time
import random

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import Consensus_Score, ProfileMatrix, RandomlyChosenKmer


def GibbsSampler(DNA, k, t, N):
    # This function finds the best scoring k-mer motif from DNA, using an iterative scheme of discarding one k-mer from a DNA at a time, and replacing it
    # with a k-mer from that DNA according to the probability computed using the Profile Matrix
    # input: DNA: a string of DNA
    #        k: integer
    #        t: integer, number of DNA strings 
    #        N: integer, number of iterations to perform Gibb's sampler
    # output: bestMotifs
    #         bestScore

    n = len(DNA[0])
    assert t == len(DNA), "number of DNA strings does not match k."

    # Generate a Motifs of k-mers by randomly picking one k-mer from each DNA
    Motifs = []
    for i in range(t):
        y = random.randint(0, n-k)
        Motifs.append(DNA[i][y:y+k])

    # score the Motifs, set it to be bestScore
    _, bestScore = Consensus_Score(Motifs)
    bestMotifs = Motifs

    # loop through the Gibb's Sampling Algorithm N times
    for j in range(N):
        random.seed(10+j)
        i = random.randint(0, t-1) # pick a random string i in DNA
        Motifs.pop(i)
        Profile = ProfileMatrix(Motifs, PseudoCounts = True)
        
        # Choose a Kmer from i-th DNA according to the probability distribution based on the Profile Matrix
        Kmer, _ = RandomlyChosenKmer(DNA[i], Profile)
        Motifs.insert(i, Kmer) # insert the randomly picked k-mer back into i-th position
        _, Score = Consensus_Score(Motifs)
        
        if Score < bestScore:
            bestMotifs = Motifs
            bestScore = Score

    return bestMotifs, bestScore
