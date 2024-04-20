import time
from datetime import timedelta
import os
import argparse
import json
import logging
from Bio import AlignIO

import matplotlib as mpl
mpl.use("Agg")


from msasa.simulated_annealing import SimulatedAnnealing
from msasa import io, plotter
from msasa.objective_functions_builder import Builder

def parse_arguments():
    parser = argparse.ArgumentParser(description='MSASA: Simulated Annealing for Multiple Sequence Alignment')
    parser.add_argument('fasta_file', type=str, help='Path to the input FASTA file', default="./data/example.fasta")
    parser.add_argument('--experiments_log_file', type=str, help='Path to the experiment JSON file', default="./data/example.fasta")
    parser.add_argument(
        "aligned_fasta_file",
        type=str,
        help="Path to the output aligned FASTA file",
        default="./data/aligned_example.fasta",
    )
    parser.add_argument('--experiments', type=int, default=1, help="Number of repetitions")
    parser.add_argument(
        "--quality_function",
        type=str,
        choices=["coincidences", "identity", "similarity_blosum62", "similarity_pam250", "similarity_gonnet", "global", "local"],
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
        "Timestamp\tIteration\tMax_Length\tTemperature\tCurrent_Score\tIteration_Score\tNew_Score\tScore_Change\tScore_Change_Perc\tBest_Score\tHistorical_Score_Improvement\tScore_Improvement_Perc\tTotal_Accepted\tTotal_Rejected\tAccepted\tNo_Changes\tAcceptance"
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

    sequence_ids, sequences = io.load_fasta(args.fasta_file)
    builder = Builder(args.match_score, args.mismatch_score, args.gap_score)
    quality_function = builder.create_quality_instance(args.quality_function)

    results_list = []
    simulated_annealing = SimulatedAnnealing(
        sequences,
        logger=None,
        extend=args.extend,
        temp=args.temperature,
        cooling_rate=args.cooling_rate,
        quality_function=quality_function,
        min_temp=args.min_temperature,
        no_changes_limit=args.max_no_changes,
        changes=args.changes,
        iteration_neighbors=args.iteration_neighbors,
    )

    for experiment_index in range(args.experiments):
        log_file = os.path.abspath(args.log_file.replace(".log", f".{experiment_index}.log"))
        print("Experiment #", experiment_index, log_file)

        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        file_handler = logging.FileHandler(log_file, mode="w")
        formatter = logging.Formatter("%(asctime)s.%(msecs)03d\t%(message)s", "%Y-%m-%d %H:%M:%S")
        file_handler.setFormatter(formatter)

        logger = logging.getLogger()
        logger.addHandler(file_handler)
        logger.setLevel(logging.INFO)

        log_header(logger)

        simulated_annealing.logger = logger

        start_time = time.time()
        aligned_sequences, final_temp, initial_score, final_score, final_iteration = simulated_annealing.anneal()
        end_time = time.time()

        elapsed_time = end_time - start_time
        readable_time = str(timedelta(seconds=elapsed_time))

        saved_alignment_file = io.save_alignment(sequence_ids, aligned_sequences, args.aligned_fasta_file, experiment_index)
        print("Output", saved_alignment_file)

        results_list.append({
            "experiment_index": experiment_index,
            "initial_temp": f"{args.temperature:.7f}",
            "final_temp": f"{final_temp:.7f}",
            "initial_score": f"{initial_score:.7f}",
            "final_score": f"{final_score:.7f}",
            "final_iteration": final_iteration,
            "output_file": saved_alignment_file,
            "elapsed_time": elapsed_time,
            "readable_time": readable_time
        })

        original_sequences = io.parse_fasta(args.fasta_file)
        aligned_sequences = io.parse_fasta(saved_alignment_file)

        assert_sequences(original_sequences, aligned_sequences)

        # Create plots
        plotter.plot_charts_from_log(log_file, plotter.plot_title(args), args.max_no_changes)
        plotter.plot_seaborn_charts_from_log(log_file, plotter.plot_title(args))

        # Read the alignment from the aligned FASTA file
        alignment = AlignIO.read(saved_alignment_file, "fasta")

        # Write the alignment to an MSF file
        AlignIO.write(alignment, saved_alignment_file + ".clustal", "clustal")

    with open(args.experiments_log_file, "w") as fp:
        json.dump(results_list, fp, indent=4)

if __name__ == "__main__":
    main()
