import os
import sys
import time
from itertools import product

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import read_data_file, SumHammingDistance
from algorithms.MedianString import MedianString


def main():
    
    data_name = "subtle_motif_dataset.txt"
    string_list = read_data_file(data_name)

    # DNA = string_list[1:]
    # k = int(string_list[0])

    DNA = string_list
    k = 8 # takes around 700 sec

    start_time = time.time()  # record start

    Median, Distance = MedianString(DNA, k)

    end_time = time.time()    # record end

    print(f"The median strings are {' '.join(Median)}, which is a distance of {Distance} from DNA.")

    print(f"Duration: {end_time - start_time:.4f} seconds")

    # data_name = "dataset_30312_1.txt"
    # string_list = read_file(data_name, dataset_dir)

    # pattern = string_list[0]
    # string = string_list[1].split()  
    # sum = SumHammingDistance(pattern, string)
    # print(sum)

if __name__ == "__main__":
    main()


