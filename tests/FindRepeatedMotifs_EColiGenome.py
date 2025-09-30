import os
import time
import sys

# Add the parent directory of this file to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from algorithms.Utils import read_data_file, FindClumps

def main():

    data_name = "E_coli.txt"
    # Read the data file, return a list of strings
    string_list = read_data_file(data_name)

    Text = string_list[0]
    # k,L,t = [int(x) for x in string_list[1].split(" ")]
    # k, L, t = [5, 500, 3]
    k, L, t = [5, 500, 3]
    # By running FindClumps on Ecoli genome, we search for short DNA words that are unusually frequent in localized regions.
    # The problem models the detection of regulatory motifs, e.g. transcription binding site, which are often clustered together near promoters.
    Unique_Patterns = FindClumps(Text, k, L, t)
    no_UP = len(Unique_Patterns)
    print(f"The {no_UP} unique {k}-mers in any {L} window with frequency >= {t} are {' '.join(Unique_Patterns)} ")


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
