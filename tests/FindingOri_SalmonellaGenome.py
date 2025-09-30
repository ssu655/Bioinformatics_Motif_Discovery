import os
import sys
import time

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import read_data_file, PlotSkew, MinimumSkew, FrequentWordsWithMismatchesAndRevComp

def main():
    
    data_name = "Salmonella_enterica.txt"
    # Read the data file, return a list of strings
    string_list = read_data_file(data_name)

    sequence = "".join(line.strip() for line in string_list if not line.startswith(">"))
    # sequence = string_list[0]

    min_index, skew = MinimumSkew(sequence)
    PlotSkew(skew, min_index)

    # the k-mers we are looking for
    k = 9
    # the number of mismatches we allow for
    d = 1
    # the half window width for replication origin
    window_width = 250
    window = sequence[min_index[0]-window_width : min_index[0]+window_width]

    PotentialDNABox, m = FrequentWordsWithMismatchesAndRevComp(window,k,d)
    print(f"The {len(PotentialDNABox)} most frequent {k}-mers, with at most {d} mismatches are {' '.join(PotentialDNABox)}. \nThey each appear {m} times.")

if __name__ == "__main__":
    # Start timer
    start_time = time.time()

    # Call main function
    main()

    # End timer
    end_time = time.time()

    # Duration in seconds
    duration = end_time - start_time
    print(f"Duration: {duration:.4f} seconds")
