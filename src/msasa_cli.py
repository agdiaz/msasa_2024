import argparse
import logging
from Bio import AlignIO

import matplotlib as mpl
mpl.use("Agg")


from msasa.simulated_annealing import SimulatedAnnealing
from msasa import io, plotter


def parse_arguments():
    parser = argparse.ArgumentParser(description='MSASA: Simulated Annealing for Multiple Sequence Alignment')
    parser.add_argument('fasta_file', type=str, help='Path to the input FASTA file', default="./data/example.fasta")
    parser.add_argument(
        "aligned_fasta_file",
        type=str,
        help="Path to the output aligned FASTA file",
        default="./data/aligned_example.fasta",
    )
    parser.add_argument(
        "--quality_function",
        type=str,
        choices=["coincidences", "identity", "similarity_blosum62", "similarity_pam250", "global", "local"],
        default="global",
        help="Quality function to evaluate alignments",
    )
    parser.add_argument(
        "--log_file",
        type=str,
        help="Path to the log file",
        default="./data/example.log",
    )
    parser.add_argument('--extend', action='store_true', help='Allow sequences to be extended beyond their original length')
    parser.add_argument('--temperature', type=float, default=5.0, help='Initial temperature for simulated annealing')
    parser.add_argument('--cooling_rate', type=float, default=0.999, help='Cooling rate for simulated annealing')
    parser.add_argument('--min_temperature', type=float, default=0.000001, help='Minimal temperature of the system')
    parser.add_argument('--max_no_changes', type=int, default=15000, help='Max consecutive no change events before early stopping')
    parser.add_argument('--match_score', type=float, default=1.0, help='')
    parser.add_argument('--mismatch_score', type=float, default=-1.0, help='')
    parser.add_argument("--gap_score", type=float, default=-10.0, help="")
    parser.add_argument("--changes", type=int, default=1, help="")
    parser.add_argument("--iteration_neighbors", type=int, default=1, help="")

    return parser.parse_args()


def log_header(logger):
    # Temporarily change the formatter to exclude timestamp
    original_formatter = logger.handlers[0].formatter
    no_timestamp_formatter = logging.Formatter("%(message)s")
    logger.handlers[0].setFormatter(no_timestamp_formatter)

    # Log the header
    logging.info(
        "Timestamp\tIteration\tMax_Length\tTemperature\tCurrent_Score\tIteration_Score\tNew_Score\tScore_Change\tScore_Change_Perc\tBest_Score\tHistorical_Score_Improvement\tScore_Improvement_Perc\tTotal_Accepted\tTotal_Rejected\tAccepted\tNo_Changes\tAcceptance\tSequences_Score_Hits\tSequences_Score_Cached\tColumn_Score_Hits"
    )

    # Revert the formatter to include timestamps for subsequent log entries
    logger.handlers[0].setFormatter(original_formatter)


def assert_sequences(original_sequences, aligned_sequences):
    # Assert that all aligned sequences have the same length
    aligned_lengths = [len(seq) for seq in aligned_sequences.values()]

    assert (len(set(aligned_lengths)) == 1), "Not all aligned sequences have the same length."
    if not len(set(aligned_lengths)) == 1:
        exit(1)

    # Assert that the length of each original sequence matches its aligned version, ignoring gaps
    for seq_id, original_seq in original_sequences.items():
        aligned_seq = aligned_sequences.get(seq_id)
        assert aligned_seq is not None, f"Aligned sequence for {seq_id} not found."
        if aligned_seq is None:
            exit(2)

        original_length = len(original_seq)
        aligned_length  = len(aligned_seq.replace("-", ""))

        assert (original_length == aligned_length), f"Length mismatch for {seq_id}: original ({original_length}) vs aligned ({aligned_length})"
        if original_length != aligned_length:
            exit(3)


def main():
    args = parse_arguments()

    logging.basicConfig(
        filename=args.log_file,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s.%(msecs)03d\t%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger = logging.getLogger()

    log_header(logger)

    sequence_ids, sequences = io.load_fasta(args.fasta_file)

    sa = SimulatedAnnealing(
        sequences,
        logger,
        extend=args.extend,
        temp=args.temperature,
        cooling_rate=args.cooling_rate,
        quality_function=args.quality_function,
        min_temp=args.min_temperature,
        no_changes_limit=args.max_no_changes,
        match_score=args.match_score,
        mismatch_score=args.mismatch_score,
        gap_penalty=args.gap_score,
        changes=args.changes,
        iteration_neighbors=args.iteration_neighbors,
    )

    sa.anneal()

    io.save_alignment(sequence_ids, sa.sequences, args.aligned_fasta_file)

    original_sequences = io.parse_fasta(args.fasta_file)
    aligned_sequences  = io.parse_fasta(args.aligned_fasta_file)

    assert_sequences(original_sequences, aligned_sequences)

    # Create plots
    plotter.plot_charts_from_log(args.log_file, plotter.plot_title(args), args.max_no_changes)

    # Read the alignment from the aligned FASTA file
    alignment = AlignIO.read(args.aligned_fasta_file, "fasta")

    # Write the alignment to an MSF file
    AlignIO.write(alignment, args.aligned_fasta_file + ".clustal", "clustal")


if __name__ == "__main__":
    main()
