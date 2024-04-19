import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def load_fasta(fasta_file):
    sequences = []
    sequence_ids = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        sequence_ids.append(record.id)

    max_length = max(len(seq) for seq in sequences)
    if np.random.rand() > 0.5:
        padded_sequences = [seq.ljust(max_length, "-") for seq in sequences]
    else:
        padded_sequences = [seq.rjust(max_length, "-") for seq in sequences]

    sequences_matrix = np.array([list([ord(c) for c in seq]) for seq in padded_sequences], dtype=np.uint8)

    return sequence_ids, sequences_matrix


def save_alignment(sequence_ids, alignment, output_file="aligned_sequences.fasta", experiment_index=1):
    seq_records = []

    for idx, seq in enumerate(alignment):
        # Convert array back to string and remove trailing gaps
        sequence_str = "".join(seq.view(dtype='S1').astype(str))
        seq_record = SeqRecord(Seq(sequence_str), id=sequence_ids[idx], description="")
        seq_records.append(seq_record)

    SeqIO.write(seq_records, output_file.replace(".aln", f".{experiment_index}.aln"), "fasta")
    return output_file.replace(".aln", f".{experiment_index}.aln")


def parse_fasta(fasta_file):
    """Parses a FASTA file and returns a dictionary of sequences."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    return sequences
