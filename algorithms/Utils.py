import os
import math
import matplotlib.pyplot as plt
from itertools import product
import random

def read_data_file(file_name):
    # read file_name in working directory
    # return a list of strings of each lines
    
    # Get the root directory of the project
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # Get the directory of the datasets folder
    dataset_dir = os.path.join(root_dir, "datasets")

    # make the full directory of the dataset file
    file_path = os.path.join(dataset_dir, file_name)

    # make an empty list for each line in dataset file
    string_list = []

    # read the dataset file line-by-line
    with open(file_path, "r") as file:  # "r" means read mode
        for line in file:
            line_str = line.strip()
            string_list.append(line_str)

    return string_list

def FrequencyTable(Text,k):
    # Generate a dictionary of the frequencies of all the k-mers in Text
    # Text is string and k is integer 
    freqMap = {}
    n = len(Text)
    for i in range(n-k):
        Pattern = Text[i:i+k]
        if Pattern not in freqMap:
            freqMap[Pattern]=1 
        else:
            freqMap[Pattern]+=1
    return freqMap

def FindClumps(Text, k, L, t):
    # this function finds the k-mer in a L-window of Text, which appears at least t times.
    Patterns = []
    n = len(Text)
    int = math.floor((n-L)/10) # compute the 1/10th interval size of genome, for reporting progress
    for i in range(n-L):
        window = Text[i:i+L]
        freqMap = FrequencyTable(window,k)
        for key in freqMap:
            if freqMap[key]>=t:
                Patterns.append(key)
        # report 10% progress
        if (i % int)==0:
            print(f"Progress {(i+1)/(n-L)*100:.1f}% of total genome.")
    # remove duplicated Patterns
    Unique_Patterns = list(set(Patterns))

    return Unique_Patterns 

def MinimumSkew(Text):
    # return the location in DNA (Text) with the minimum skew
    # the minimum skew index is indicative of where the origin of replication is in a genome
    # skew(i+1) = no. of G - no. of C from Text[0:i], or the difference between accumulated G - C 
    # Outputs: min_index is a list of indices (int) which corresponds to min_skew
    #          skew is a list of (G-C) numbers in integers
    n = len(Text)
    skew = [0]
    for i in range(0,n):
        if Text[i]=="C":
            skew.append(skew[i]-1)
        elif Text[i]=="G":
            skew.append(skew[i]+1)
        else:
            skew.append(skew[i])
    min_skew = min(skew)
    min_index = [i for i, val in enumerate(skew) if val == min_skew]
    return min_index, skew

def PlotSkew(skew, min_index):
    # Create the plot
    plt.plot(skew)
    plt.plot(min_index, [skew[i] for i in min_index], marker='o')

    # Add labels
    plt.xlabel("Distance along genome")
    plt.ylabel("#G - #C")

    # Add a title (optional)
    plt.title("Plot the skew (#G-#C) along the genome")

    # Show the plot
    plt.show()

    return

def ReverseComplement(Text):
    # return the reverse complement of a DNA strand (Text)
    # A-T, C-G
    RC = ''
    n = len(Text)
    for i in range(n-1,-1,-1):
        if Text[i]=="A":
            RC += "T"
        elif Text[i]=="T":
            RC += "A"
        elif Text[i]=="C":
            RC += "G"
        elif Text[i]=="G":
            RC += "C"
        else: # wildcard nucleotide
            RC += Text[i]

    return RC

def ImmediateNeighbors(Pattern):
    # Generate the 1-neighborhood of Pattern (a k-mer), the set of all k-mers who are the same or 1-off Pattern
    # or Hamming distance no more than 1
    # output Neighborhood is a list of unique strings
    Neighborhood = [Pattern,]
    nucleotides = ["A", "T", "C", "G"]
    for i in range(len(Pattern)):
        for nt in nucleotides:
            if nt != Pattern[i]:
                neighbor = Pattern[:i] + nt + Pattern[i+1:]
                Neighborhood.append(neighbor)
    return Neighborhood

def Neighbors(Pattern, d):
    # Generate the d-neighborhood of Pattern (a k-mer), the set of all k-mers who are less than d off Pattern
    # or Hamming distance no more than d
    # output Neighborhood is a list of unique strings
    Neighborhood = []
    i_Neighbors = [Pattern,]
    
    # just need to repeatedly compute the immediate neighborhoods of a pattern, n times, to get the n-neighbors
    for i in range(d):
        for n in i_Neighbors:
            Neighbors = ImmediateNeighbors(n)
            Neighborhood = Neighborhood + Neighbors
        # there are likely duplicates between each immediate neighborhoods, so make it unique
        i_Neighbors = list(set(Neighborhood))
    Neighborhood = i_Neighbors    

    return Neighborhood

def FrequentWordsWithMismatchesAndRevComp(Text,k,d):
    # The function takes in a DNA (Text) and compute the most frequent k-mers, with at most d mismatch nucleotides.
    # outputs a list of most frequent k-mers
    Patterns = []
    freqMap = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        Pattern_rc = ReverseComplement(Pattern)
        
        # make a huge joint list of neighborhoods from the pattern and Rev Comp of pattern
        # hopefully there are no overlaps. You'd think for long enough patterns, there shouldn't be!
        # should probably make the overall neighborhood unique actually
        Neighborhood = Neighbors(Pattern,d) + Neighbors(Pattern_rc,d)

        for j in range(len(Neighborhood)):
            # loop through each pattern in neighborhoods
            # the idea is, the pattern is the k-neighbor of each of its k-neighbors
            # so each time a pattern shows up, all of it's k-neighbors increase occurance count by 1
            neighbor = Neighborhood[j]
            if neighbor not in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] += 1
    m = max(freqMap.values())
    for key in freqMap:
        if freqMap[key] == m:
            Patterns.append(key)
    
    return Patterns, m

def Consensus_Score(Motifs):
    # This function computes the Consensus AND the Score 
    # consensus: the k-mer consisting of the most popular nucleotides in a list of Motifs, 
    # score : the number of unpopular nucleotides in all of Motifs, wildcard nucleotides ignored
    # 
    # input: Motifs: a list of strings of same length
    # output: consensus: a string
    #         Score: an integer number

    consensus = "" 
    score = 0
    bases = ["A","T","C","G"]
    
    n = len(Motifs) # length of each motif
    k = len(Motifs[0])

    for i in range(k):
        # Count nucleotides at position i 
        counts = {}
        for Motif in Motifs:
            nucleotide = Motif[i]
            if nucleotide in counts:
                counts[nucleotide] += 1
            else:
                counts[nucleotide] = 1
        # Find the nucleotide with the maximum count
        max_count = 0
        max_nuc = ''
        for nuc, count in counts.items():
            if count > max_count:
                max_count = count
                max_nuc = nuc
        consensus += max_nuc
        # count the total number of bases that are not max_nuc, do not count any wildcard nucleotides
        non_pop_nuc_count = 0
        for j in range(4):
            if bases[j] != max_nuc and bases[j] in counts:
                non_pop_nuc_count += counts[bases[j]]
        # adding non popular nucleotide count to the score
        score += non_pop_nuc_count

    return consensus, score

def ProfileMatrix(Motifs, PseudoCounts = False):
    # This function returns the Profile Matrix of Motifs
    # argument: PseudoCounts: Boolean variable, determining whether to apply Laplace's Rule of Succession, substituting zeros with small numbers called pseudocounts
    # input: Motifs: a list of k-mer strings
    # output: Profile: a dict list of float, of dimention [4, k]
    
    profile_keys = ["A", "C", "G", "T"]
    profile = {}
    for i in profile_keys:
        profile[i] = [] 
    t = len(Motifs)
    k = len(Motifs[0])
    for i in range(k):
        column = [motif[i] for motif in Motifs]
        for j in range(len(profile_keys)):
            count = column.count(profile_keys[j])
            profile[profile_keys[j]].append(count)
    if PseudoCounts:
        t += 4 
        for i in profile_keys:
            profile[i] = [j+1 for j in profile[i]]
    for i in profile_keys:
        profile[i] = [j/t for j in profile[i]]

    return profile

def MostProbableKmer(Text, Profile, k ,IgnoreWildCard = True):
    # This function computes the most probable k-mer among all k-mers from Text according to the Profile matrix
    # Input: Text: a string of length at least k
    #        Profile: a dict of 4 list of float numbers, list length k, whose elements are less than 1 and each column sums to 1
    #        k : an integer
    # Output: highest_kmer: a string of most probably kmer
    #         highest_prob: a float number of probability

    assert k == len(Profile["A"]), "number of columns in Profile does not equal k."

    profile_keys = ["A", "C", "G", "T"]
    n = len(Text)
    
    # initiate the highest probability number, and most probable kmer string
    highest_prob = 0
    highest_kmer = Text[0:k]

    for i in range(k,n): # loop through all k-mers in Text
        kmer = Text[i-k:i]

        # If IgnoreWildCard is true, we skip any kmers with wildcards, as they are ambiguous
        if IgnoreWildCard:
            if '*' in kmer or 'N' in kmer:
                continue

        # set the initial probability before considering any nucleotides in kmer
        prob = 1 

        for j in range(k): # loop through the k-mer
            nt = kmer[j]

            # the j element entry of Profile[nt]
            if nt in profile_keys:
                prob *= Profile[nt][j] 
            else: # if nt is a wildcard, count it as any of the 4 nucleotides
                prob *= 1

            # if the current probability, at position j, is already lower than the most probably kmer so far, just give up this kmer
            # if (prob < highest_prob) and i!=0: 
            #   break
        # update the highest probability to new prob, and update the best kmer
        if prob > highest_prob:
            highest_prob = prob
            highest_kmer = kmer

    return highest_kmer, highest_prob

def HammingDistance(Text1, Text2):
    # Hamming distance is the number of mismatching characters between two strings of same length

    # check the lengths of two strings match
    assert len(Text1) == len(Text2), f"Lengths differ: {len(Text1)} vs {len(Text2)}"

    n = len(Text1)
    count = 0
    for i in range(n):
        if Text1[i]!=Text2[i]:
            count += 1

    return count

def MinHammingDistance(Pattern, DNA):
    # This function computes the minimum Hamming Distance between a k-mer (Pattern) and a DNA segment (DNA)
    # inputs: Pattern: the k-mer string of nucleotides
    #         DNA: the DNA segment that's same or longer than len(Pattern)
    # output: d: minimum Hamming Distance, integer

    n = len(DNA)
    k = len(Pattern)
    min_d = n # set minimum distance to be the largest possible value
    for i in range(n-k):
        d = HammingDistance(DNA[i:i+k],Pattern)
        if d < min_d:
            min_d = d
    return min_d

def SumHammingDistance(Pattern, DNA):
    # This function computes the sum of minimum Hamming Distances between a k-mer (Pattern) and a list of DNA_segment
    # inputs: Pattern: the k-mer string of nucleotides
    #         DNA: a list of DNA segment that's same or longer than len(Pattern)
    # output: sum: sum of minimum Hamming Distance, integer
    
    sum = 0
    for i in range(len(DNA)):
        sum += MinHammingDistance(Pattern, DNA[i])
    return sum

def K_merSpace(k):
    # This function returns a list of all possible k-mers from 'AAA..A' to 'TTT..T'
    # input: k: integer
    # output: k_mers: a list of strings

    nucleotides = ["A", "T", "C", "G"]
    kmers = [''.join(p) for p in product(nucleotides, repeat=k)]
    return kmers

def MostProbableMotifs(DNA, Profile):
    # This function returns the most probable Motifs in DNA given the profile Matrix Profile.
    # this function mostly uses the MostProbableKmer function
    # inputs: DNA: list of strings
    #         Profile: a dict of lists
    # output: Motifs: a list of k-mers

    t = len(DNA)
    k = len(Profile['A'])
    Motifs = []
    for i in range(t):
        Kmer, _ = MostProbableKmer(DNA[i], Profile, k)
        Motifs.append(Kmer)

    return Motifs

def Weighted_Choice(weights):
    # This function returns a integer from 0 to len(weights)-1
    # input: weights: a list of floats that sum to one
    # output: an integer between 0 and len(weights)-1 inclusive
    r = random.random()  # uniform [0,1)
    cumulative = 0.0
    for i, w in enumerate(weights):
        cumulative += w
        if r < cumulative:
            return i

def RandomlyChosenKmer(Text, Profile):
    # This function computes a randomly chosen k-mer from all k-mers from Text, according to the probability distribution computed from Profile matrix
    # Input: Text: a string of length at least k
    #        Profile: a dict of 4 list of float numbers, list length k, whose elements are less than 1 and each column sums to 1
    # Output: kmer
    #         prob: a float number of probability

    n = len(Text)
    k = len(Profile["A"])

    # initiate the highest probability number, and most probable kmer string
    prob_list = []

    for i in range(k,n): # loop through all k-mers in Text
        kmer = Text[i-k:i]
        prob = 1 # beginning with 1, compute the probability of this kmer

        for j in range(k): # loop through the k-mer
            nt = kmer[j]

            # the j element entry of Profile[nt]
            prob *= Profile[nt][j] 

        prob_list.append(prob)
        
    total_prob = sum(prob_list)
    prob_list2 = [i/total_prob for i in prob_list]

    # choose an index based on the probability distributed computed above
    i = Weighted_Choice(prob_list2)
    
    kmer = Text[i:i+k]
    prob = prob_list[i]

    return kmer, prob