# PROBLEM 1: Read a FASTA file, and output any region between a start (‘ATG’) a>
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def find_open_reading_frame(sequence):
    rframe = [] # List to store found ORFs
    seq_len = len(sequence)

    # Search for ORFs in the sequence
    for frame in range(3):
        for start in range(frame, seq_len, 3):
            codon = sequence[start:start+3]
            if codon == 'ATG': # Start codon found
                for end in range(start + 3, seq_len, 3):
                    stop_codon = sequence[end:end+3]
                    if stop_codon in ['TAA', 'TAG', 'TGA']:
                        orf  = sequence[start:end+3] # Get the ORF
                        rframe.append((start, end+3, orf))  # Store start, end,>
                        break
    return rframe

def main(input_file):
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)
        rframe = find_open_reading_frame(sequence)
        
       	for start, end, orf in rframe:
            print(f">ORF_{start}_{end}")
            print(orf)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find ORFs in a FASTA file")
    parser.add_argument("input_file", help="Input FASTA file")
    args = parser.parse_args()
    
    main(args.input_file)
