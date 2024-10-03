import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

def translate_orfs(orfs):
    """Translates a list of ORFs to protein sequences."""
    protein_strings = set()
    for _, _, _, orf in tqdm(orfs, desc="Translating ORFs"):  # Ajustado a cuatro variables
        protein = Seq(orf).translate(to_stop=True)
        protein_strings.add(str(protein))
    return protein_strings

def find_ORF(sequence):
    """Finds all ORFs in the provided DNA sequence."""
    ORFs = set()
    n = len(sequence)
    for strand, seq in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            for start in range(frame, n, 3):
                if seq[start:start + 3] == 'ATG':
                    for end in range(start + 3, n, 3):
                        if seq[end:end + 3] in ['TAA', 'TAG', 'TGA']:
                            ORF = seq[start:end + 3]
                            ORFs.add((start, end + 3, strand, ORF))
                            break
    return ORFs

def main(input_file, output_file):
    """Main function to process the input FASTA file and output protein sequences."""
    with open(output_file, "a") as f:
        for record in SeqIO.parse(input_file, "fasta"):
            sequence = Seq(str(record.seq))  # Convert to Seq object
            orfs = find_ORF(sequence)

            protein_strings = translate_orfs(orfs)

            # Save protein sequences to the output file
            for protein in protein_strings:
                f.write(f"{protein}\n")

            # Optionally print the ORFs to the console
            for start, end, _, orf in orfs:
                print(f">ORF_{start}_{end}")
                print(orf)

        print(f'Job for file {input_file} finished')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find ORFs in a FASTA file")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("output_file", help="Output file for protein sequences")
    args = parser.parse_args()

    main(args.input_file, args.output_file)
