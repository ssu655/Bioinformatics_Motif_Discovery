import os
import sys
from collections import Counter

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.GreedyMotifSearch import GreedyMotifSearch
from algorithms.Utils import read_data_file

def main():
    # Subtle motif dataset contains wildcard *, 
    # special treatment is implemented in Consesus_Score, and MostProbableKmer methods. We assume * can be any nucleotides there.
    # but it is possible to ignore any k-mers with * completely.
    data_name = "subtle_motif_dataset.txt" 
    # Read the data file
    string_list = read_data_file(data_name)

    # k, d = [int(i) for i in string_list[0].split()]
    k = 20
    DNA = string_list[1].split()
    DNA = string_list

    motif, score = GreedyMotifSearch(DNA, k, PseudoCounts = True, IgnoreWildCard = False) 
    # Note The PseudoCounts = True version of the algorithm somehow never past the test. Something is wrong

    print(f"The best Motifs are {' '.join(motif)} with a score of {score}")

if __name__ == "__main__":
    main()