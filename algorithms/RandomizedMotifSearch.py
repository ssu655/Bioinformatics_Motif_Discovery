import os
import sys
import time
import random


# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import ProfileMatrix, Consensus_Score, MostProbableMotifs


def RandomizedMotifSearch(DNA, k, t):
    # This function generates a randomly selected motif from each string in DNA, and forms a profile matrix, which is used to find a motif, so on and so forth,
    # eventually hoping to land on a final best, lowest scored motif 
    # inputs: DNA: list of strings
    #         k: integer
    #         t: integer, like time, number of different randomized initial Motifs
    # output: Motifs
    #         bestScore: integer

    n = len(DNA[0])
    d = len(DNA)
    
    bestMotifs = []
    verybestScore = n*k

    for j in range(t):
        bestScore = n*k
        Motifs = []
        for i in range(d):
            random.seed(1+i*2+t)
            y = random.randint(0, n-k)
            Motifs.append(DNA[i][y:y+k])
        _, Score = Consensus_Score(Motifs)

        while Score < bestScore:
            bestScore = Score
            Profile = ProfileMatrix(Motifs, PseudoCounts = True)
            Motifs = MostProbableMotifs(DNA, Profile)
            _, Score = Consensus_Score(Motifs)
        if bestScore < verybestScore:
            verybestScore = bestScore
            bestMotifs = Motifs
    return bestMotifs, verybestScore

