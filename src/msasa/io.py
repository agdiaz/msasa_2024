import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import matplotlib as mpl

mpl.use("Agg")

def load_fasta(fasta_file):
    sequences = []
    sequence_ids = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        sequence_ids.append(record.id)

    max_length = max(len(seq) for seq in sequences)
    padded_sequences = [seq.ljust(max_length, "-") for seq in sequences]

    sequences_matrix = np.array([list(seq) for seq in padded_sequences], dtype='S1')

    return sequence_ids, sequences_matrix


def save_alignment(sequence_ids, alignment, output_file="aligned_sequences.fasta"):
    seq_records = []

    for idx, seq in enumerate(alignment):
        # Convert array back to string and remove trailing gaps
        sequence_str = "".join(seq.astype(str))
        seq_record = SeqRecord(Seq(sequence_str), id=sequence_ids[idx], description="")
        seq_records.append(seq_record)

    SeqIO.write(seq_records, output_file, "fasta")


def parse_fasta(fasta_file):
    """Parses a FASTA file and returns a dictionary of sequences."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    return sequences
