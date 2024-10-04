# Week-3-4 (some of the problems were solved with the help of some classmates)

LLM used: CHAT GPT-4o mini
## Step 1: Initialize BioPython

ssh alboeis@ilogin.ibex.kaust.edu.sa
pip install biopython
## Step 2: Create Git Repository
```bash 
mkdir Week34
cd Week34/
git init
nano gene_finder.py
```
## Problem 1: Finding reading frames
```bash
python gene_finder.py test1.fasta > output1.txt
git git add gene_finder.py output1.txt
git commit -m "Problem 1"
```
gene_finder.py and output1.txt

## Problem 2: Finding reading frames including the reverse complements

```bash
touch gene_finder_RvCom.py
nano gene_finder_RvCom.py 
python gene_finder_RvCom.py test1.fasta > output2.txt
git add gene_finder_RvCom.py output2.txt
git commit -m "Solving Problem 2"
```
gene_finder_RvCom.py and output2.txt

## Problem 3: Solve the Open Reading Frame problem on Rosalind
```bash
nano Rosalind_orf.fna #Here we add the data given by Rosalind example
nano Rosalind_problem72.py
python3 Rosalind_problem72.py Rosalind_orf.fna
nano C_Output_Proteins.txt
```
Rosalind_problem72.py

## Problem 4: Finding reading frames including the reverse complements and applying code to all 14 downloaded genomes
```bash
nano gene_finder_adapted.py 
find /home/alboeis/ncbi_dataset/data -type f -name "*GCF*.fna" -exec python gene_finder_adapted.py {} all_orfs_adapted.txt \;
git add gene_finder_adapted.py  all_orfs_AminoAcids.txt
git commit -m "Saving files of problem 4"
```
gene_finder_adapted.py all_orfs_AminoAcids.txt

NUMBER OF LINES: 765037 all_orfs_AminoAcids.txt

## Problem 5: Implementing gene finder with length filter

```bash
nano gene_filter_new.py
find /home/alboeis/ncbi_dataset/data -type f -name "*GCF*.fna" -exec python gene_filter_new.py {} output_filter_new.txt -l 100 \;
git add gene_filter_new.py output_filter_new.txt
git commit -m "Saving files of problem 5"
```
gene_filter_new.py output_filter_new.txt

NUMBER OF LINES: 162136 output_filter_new.txt 

## Problem 6: Implementing gene finder with length, Ribosome binding site (rbs) and rbs type filter

```bash
nano gene_RBS.py
find /home/alboeis/ncbi_dataset/data -type f -name "*GCF*.fna" -exec python gene_RBS.py {} output_RBS.txt -l 100 -r AGGAGG -u 20 \;
git add gene_RBS.py gene_RBS.py
git commit -m "Files of Problem 6"
git push -u origin main
```
gene_RBS.py and output_RBS.txt

NUMBER OF LINES: 900 output_RBS.txt (7.9 K)
