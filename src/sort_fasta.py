import sys
from Bio import SeqIO

def sort_fasta(fasta_file, ref_file, output_file):
    # Read the reference order
    order_records = list(SeqIO.parse(ref_file, "fasta"))
    order = [record.id for record in order_records]

    # Load the sequences from the FASTA file
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Sort the sequences according to the reference file
    sorted_records = [record_dict[id] for id in order if id in record_dict]

    # Write sorted records to an output file
    SeqIO.write(sorted_records, output_file, "fasta")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python sort_fasta.py <fasta_file> <ref_file> <output_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    ref_file = sys.argv[2]
    output_file = sys.argv[3]

    sort_fasta(fasta_file, ref_file, output_file)
