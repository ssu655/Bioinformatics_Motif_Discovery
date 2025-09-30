import os
import sys

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import SumHammingDistance, K_merSpace

def MedianString(DNA, k):
    # This function is a brute force solution to the Median String problem, which finds the k-mers with minimum distance to the set of DNA segment
    # input: DNA: set of DNA string segments
    #        k: integer
    # output: median: a list of k-mer strings
    #         distance: an integer, the minimum distance between the median strings and DNA
    #         note that there might exist more than one strings with minimum distance, this algorithm returns all such string encountered
    #         This is related to our measure of distance, if we use entropy as a minimization measure, we will get the one and only median string
    #         because entropy is decimal
    t = len(DNA) # num of strings
    n = len(DNA[0]) # length of DNA segment

    median = []
    distance = n*t # Maximum distance possible
    for pattern in K_merSpace(k):
        d = SumHammingDistance(pattern, DNA)
        if distance > d:
            median = []
            median.append(pattern)
            distance = d
        elif distance == d:
            median.append(pattern)
    return median, distance
    