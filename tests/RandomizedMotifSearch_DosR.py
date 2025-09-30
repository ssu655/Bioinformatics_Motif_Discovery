import os
import sys
import time

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import read_data_file
from algorithms.RandomizedMotifSearch import RandomizedMotifSearch

def main():
    data_name = "DosR.txt"
    # Read the data file
    string_list = read_data_file(data_name)

    # DNA = string_list[1].split()
    # k, d = [int(i) for i in string_list[0].split(" ")]
    t = 1000 
    k = 20
    DNA = string_list

    start_time = time.time()  # record start
    Motifs,bestScore = RandomizedMotifSearch(DNA, k, t)
    end_time = time.time()    # record end

    print(f"The best motifs, are {' '.join(Motifs)}, the best Score is {bestScore}")
    print(f"Duration: {end_time - start_time:.4f} seconds")

if __name__ == "__main__":
    main()


