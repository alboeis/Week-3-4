import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def find_ORFs(sequence):
    ORFs = []
    n = len(sequence)

    for strand, seq in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            for start in range(frame, n, 3):
                if seq[start:start+3] == 'ATG':  # Start codon
                    for end in range(start + 3, n, 3):
                        if seq[end:end+3] in ['TAA', 'TAG', 'TGA']:  # Stop codons
                            ORF = seq[start:end+3]
                            ORFs.append(ORF)  # Only store the ORF
                            break
    return ORFs

def main(input_file):
    proteins = []
    unique_proteins = set()

    # Parse the input FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = record.seq
        ORFs = find_ORFs(sequence)

        for idx, orf in enumerate(ORFs):
            protein = orf.translate()  # Translate the ORF to protein
            if protein:  # Only add non-empty proteins
                protein_str = str(protein).rstrip('*')  # Remove trailing '*'
                proteins.append(f">Protein_{idx + 1}\n{protein_str}")
                unique_proteins.add(protein_str)  # Use the version without '*'

    # Save the first file with headers
    with open('C_Output_Proteins_with_headers.txt', 'w') as f:
        f.write('\n'.join(proteins))

    # Save the second file with proteins on separate lines
    with open('C_Output_Proteins.txt', 'w') as f:
        for protein in unique_proteins:
            f.write(f"{protein}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find ORFs in a FASTA file")
    parser.add_argument("input_file", help="Input FASTA file")
    args = parser.parse_args()

    main(args.input_file)
