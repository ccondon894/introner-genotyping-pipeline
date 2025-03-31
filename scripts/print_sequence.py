import argparse
from pyfaidx import Fasta

def extract_sequence(fasta_file, chrom, start, end, output_file=None):
    """
    Extracts a sequence from a FASTA file based on chromosome, start, and end positions.

    :param fasta_file: Path to the input FASTA file.
    :param chrom: Chromosome name.
    :param start: Start position (1-based).
    :param end: End position (1-based, inclusive).
    :param output_file: Path to the output file (optional).
    """
    # Load the FASTA file
    fasta = Fasta(fasta_file)

    try:
        # Extract the sequence
        sequence = fasta[chrom][start - 1:end].seq  # Convert to 0-based indexing

        # Create the FASTA header
        header = f">{chrom}:{start}-{end}"

        # Print or write the sequence
        if output_file:
            with open(output_file, "w") as f:
                f.write(f"{header}\n{sequence}\n")
            print(f"Sequence written to {output_file}")
        else:
            print(f"{header}\n{sequence}")
    except KeyError:
        print(f"Chromosome '{chrom}' not found in the FASTA file.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Extract a sequence from a FASTA file based on chrom, start, and end.")
    parser.add_argument("fasta", help="Path to the FASTA file")
    parser.add_argument("chrom", help="Chromosome name")
    parser.add_argument("start", type=int, help="Start position (1-based)")
    parser.add_argument("end", type=int, help="End position (1-based, inclusive)")
    parser.add_argument("-o", "--output", help="Output file to save the extracted sequence (optional)")

    args = parser.parse_args()

    # Call the extraction function
    extract_sequence(args.fasta, args.chrom, args.start, args.end, args.output)