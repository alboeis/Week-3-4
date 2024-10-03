# Week-3-4
WEEK 3 & 4

LLM used: CHAT GPT-4o mini
## Step 1: Initialize BioPython

ssh alboeis@ilogin.ibex.kaust.edu.sa
pip install biopython
## Step 2: Create Git Repository
```bash 
mkdir Week4
cd Week4/
git init
touch gene_finder.py README.md
nano gene_finder.py
```
## Problem 1: Finding reading frames

Read a FASTA file, and output any region between a start (‘ATG’) and stop codon (‘TAA’, ‘TAG’, ‘TGA’) in that FASTA file

We need to run the code of "gene_finder.py" to find the open reading frames and I used the file "test1.fasta" to test the code.
```bash
python gene_finder.py test1.fasta > output1.txt
git git add gene_finder.py output1.txt
git commit -m "Solve problem 1"
```

SUBMITTED FILES: gene_finder.py and output1.txt


## Problem 2: Finding reading frames including the reverse complements

We need to run the code of "gene_finder_RvCom.py" to find the open reading frames and I used the file "test1.fasta" to test the code.
```bash
touch gene_finder_RvCom.py
nano gene_finder_RvCom.py 
python gene_finder_RvCom.py test1.fasta > output2.txt
git add gene_finder_RvCom.py output2.txt
git commit -m "Solving Problem 2"
```
SUBMITTED FILES: gene_finder_RvCom.py and output2.txt

## Problem 3: Solve the Open Reading Frame problem on Rosalind (Problem 72)
```bash
nano Rosalind_orf.fna #Here we add the data given by Rosalind example
nano Rosalind_problem72.py
python3 Rosalind_problem72.py Rosalind_orf.fna
nano C_Output_Proteins.txt
```
SUBMITTED FILES: Rosalind_problem72.py

## Problem 4: Finding reading frames including the reverse complements and applying code to all 14 downloaded genomes
```bash
nano gene_finder_adapted.py 
find /home/alboeis/ncbi_dataset/data -type f -name "*GCF*.fna" -exec python gene_finder_adapted.py {} all_orfs_adapted.txt \;
git add gene_finder_adapted.py  all_orfs_AminoAcids.txt
git commit -m "Saving files of problem 4"
```
FILES SUBMITTED: gene_finder_adapted.py all_orfs_AminoAcids.txt (A file containing all the ORF for the 14 downloaded genomes but transduced to aminoacids)

NUMBER OF LINES: 765037 all_orfs_AminoAcids.txt (60M)

## Problem 5: Implementing gene finder with length filter

We need to run the code of "GeneFinderFilter.py" to implement a filter by length that discards short ORFs that are unlikely to be functional genes (e.g., less than 100 codon)
```bash
nano gene_filter_new.py
find /home/alboeis/ncbi_dataset/data -type f -name "*GCF*.fna" -exec python gene_filter_new.py {} output_filter_new.txt -l 100 \;
git add gene_filter_new.py output_filter_new.txt
git commit -m "Saving files of problem 5"
```
FILES SUBMITTED: gene_filter_new.py output_filter_new.txt

NUMBER OF LINES: 162136 output_filter_new.txt (45 M)

## Problem 6: Implementing gene finder with length, Ribosome binding site (rbs) and rbs type filter

We need to run the code of "GeneFinderRBS.py" to filter all predicted ORFs based on whether they contain a Shine-Dalgarno sequence (AGGAGG) up to 20bp upstream of the start codon.
```bash
nano gene_RBS.py
find /home/alboeis/ncbi_dataset/data -type f -name "*GCF*.fna" -exec python gene_RBS.py {} output_RBS.txt -l 100 -r AGGAGG -u 20 \;
git add gene_RBS.py gene_RBS.py
git commit -m "Files of Problem 6"
git push -u origin main
```
FILES SUBMITTED: gene_RBS.py and output_RBS.txt

NUMBER OF LINES: 900 output_RBS.txt (7.9 K)
