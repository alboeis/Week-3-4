import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def find_rframes(sequence):
    rframes = []  # List to store reading frames
    seq_len = len(sequence)

    # Check both strands. Loop over both the forward and reverse strands of the DNA sequence.
    for strand, seq in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            for start in range(frame, seq_len, 3):
                if seq[start:start+3] == 'ATG':  # Start codon
                    for end in range(start + 3, seq_len, 3):
                        if seq[end:end+3] in ['TAA', 'TAG', 'TGA']:  # Stop codons
                            rframe = seq[start:end+3]  # Get rframe
                            rframes.append((start, end+3, strand, rframe))
                            break
    return rframes


def main(input_file):
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = record.seq  # Get the sequence
        rframes = find_rframes(sequence)  # Find rframes

        # Print rframes with positions and strand
        for start, end, strand, rframe in rframes:
            strand_symbol = '+' if strand == 1 else '-'
            print(f">RFRAME_{start}_{end}_{strand_symbol}")
            print(rframe)

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Find reading frames in a FASTA file")
    parser.add_argument("input_file", help="Input FASTA file")
    args = parser.parse_args()

    main(args.input_file)
