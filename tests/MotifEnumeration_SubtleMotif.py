import os
import sys
import time

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import read_data_file
from algorithms.MotifEnumeration import MotifEnumeration

def main():
    
    data_name = "subtle_motif_dataset.txt"
    # This is an artificial dataset, constructed to be challenging for motif dicovery algorithm. and very sensitive to algorithm logics
    # Read the data file
    string_list = read_data_file(data_name)
    
    # DNA = string_list[1:]
    # k = int(string_list[0])
    # d = 1

    DNA = string_list
    k, d = [10,2] #[int(i) for i in string_list[0].split(" ")]

    start_time = time.time()  # record start

    Patterns = MotifEnumeration(DNA, k, d)

    end_time = time.time()    # record end

    print(f"The {len(Patterns)} ({k},{d})-motifs, are {' '.join(Patterns)}")
    print(f"Duration: {end_time - start_time:.4f} seconds")

if __name__ == "__main__":
    main()


