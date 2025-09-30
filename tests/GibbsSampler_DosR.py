import os
import sys
import time

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import read_data_file
from algorithms.GibbsSampler import GibbsSampler

def main():

    data_name = "DosR.txt"
    # Read the data file
    string_list = read_data_file(data_name)

    # DNA = string_list[1].split()
    # k, t, N = [int(i) for i in string_list[0].split(" ")]
    DNA = string_list
    k = 20
    t = 10
    N = 1000

    start_time = time.time()  # record start
    bestScore = 100
    bestMotifs = []
    for i in range(20):
        Motifs,Score = GibbsSampler(DNA, k, t, N)
        if Score < bestScore:
            bestMotifs = Motifs
            bestScore = Score
    end_time = time.time()    # record end
    n = len(bestMotifs)

    print(f"The {n} best motifs, are {' '.join(bestMotifs)}, the best Score is {bestScore}")
    print(f"Duration: {end_time - start_time:.4f} seconds")

if __name__ == "__main__":
    main()


