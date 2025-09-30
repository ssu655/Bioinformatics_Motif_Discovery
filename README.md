# Bioinformatics_Motif_Discovery
This repository contains my implementations of some classic motif discovery algorithms. 

The goal is to showcase my ability to process genomic data and implement computational algorithms for finding motifs in DNA sequences.

The contents are adapted from the bioinformatics MOOC “Finding Hidden Messages in DNA” 

————

Project Structure

Bioinformatics-motif-discovery/
- algorithms/ 
	- \_\_init\_\_.py  
	- Utils.py 
	- GibbsSampler.py 
	- GreedyMotifSearch.py 
 	- MedianString.py 
	- MotifEnumeration.py 
	- RandomizedMotifSearch.py 
- datasets/
	- DosR.txt
	- E_coli.txt
	- Salmonella_enterica.txt
	- subtle_motif_dataset.txt
- tests/
	- FindRepeatedMotifs_EColiGenome.py
	- FindingOri_SalmonellaGenome.py
	- GreedyMotifSearch_SubtleMotif.py
	- MotifEnumeration_SubtleMotif.py
	- MedianString_SubtleMotif.py
	- RandomizedMotifSearch_DosR.py
	- GibbsSampler_DosR.py
- notebooks/
	- Demo.ipynb
- README.md
- .gitignore
- requirements.txt
- LICENSE

————

Algorithms Implemented

- Clump Finding  		– detecting overrepresented k-mers in genomic windows  
- Replication Origin Finder  – discovering the region of Replication Origin  
- Motif Enumeration 		– exhaustive motif finding 
- Median String 		– brute force optimal motif search  
- Greedy Motif Search 	– heuristic motif search  
- Randomized Motif Search – stochastic search for motifs  
- Gibbs Sampler 		– probabilistic motif discovery  

————

Example Datasets of Implementations 

- E. coli genome – find repeated motifs in clumps  
- Salmonella genome – locate the origin of replication  
- DosR regulon (M. tuberculosis)  – discover regulatory motifs using Gibbs Sampling and Randomized Motif Search  
- Subtle motif datasete - an artifically constructed DNA dataset with a subtle hidden motif

A full demonstration is provided in notebooks/Demo.ipynb

————

Installation

Clone the repository:
```
git clone https://github.com/ssu655/bioinformatics-motif-discovery.git
cd bioinformatics-motif-discovery
```

Create and activate the virtual environment
```
python -m venv venv
source venv/bin/activate   # On Mac/Linux
venv\Scripts\activate      # On Windows
pip install -r requirements.txt
```

Run the demo jupyter notebook
```
jupyter notebook notebooks/Demo.ipynb
```

————

Dependencies

Python 3.8+

Jupyter Notebook

See requirements.txt for the full list.
