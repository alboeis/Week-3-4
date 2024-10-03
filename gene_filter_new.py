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

def main(input_file, output_file, min_length):
    """Main function to process the input FASTA file and output protein sequences."""
    min_length_codons = min_length * 3  # Convert codons to nucleotides

    with open(output_file, "a") as f:
        for record in SeqIO.parse(input_file, "fasta"):
            sequence = Seq(str(record.seq))  # Convert to Seq object
            orfs = find_ORF(sequence)

            # Filter out ORFs shorter than the specified minimum length
            filtered_orfs = {orf for orf in orfs if len(orf[3]) >= min_length_codons}

            if not filtered_orfs:  # If no ORFs pass the filter, continue to the next record
                continue

            print(f">ORF")
            print(filtered_orfs)

            protein_strings = translate_orfs(filtered_orfs)

            # Save protein sequences to the output file
            for protein in protein_strings:
                f.write(f"{protein}\n")

            # Optionally print the ORFs to the console
            for start, end, _, orf in filtered_orfs:
                print(f">ORF_{start}_{end}")
                print(orf)

        print(f'Job for file {input_file} finished')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find ORFs in a FASTA file")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("output_file", help="Output file for protein sequences")
    parser.add_argument("-l", "--min_length", type=int, default=100,
                        help="Minimum ORF length in codons (default: 100)")
    args = parser.parse_args()

    main(args.input_file, args.output_file, args.min_length)
