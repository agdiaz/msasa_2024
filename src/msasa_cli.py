import os
import numpy as np
import argparse
from Bio import SeqIO
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import pandas as pd
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter

from msasa.sa import SimulatedAnnealing


def ms_formatter(a, b):
    t = mdates.num2date(a)
    ms = str(t.microsecond)[:3]
    res = f"{t.hour:02}:{t.minute:02}:{t.second:02}.{ms}"
    return res


def load_fasta(fasta_file):
    """
    Loads sequences from a FASTA file and returns them in a NumPy matrix, ensuring all sequences are of the same length.

    Parameters:
    - fasta_file (str): The path to the input FASTA file.

    Returns:
    - tuple: (sequence_ids, sequences_matrix) where sequence_ids is a list of sequence IDs and sequences_matrix is a NumPy array.
    """

    sequences = []
    sequence_ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        sequence_ids.append(record.id)

    max_length = max(len(seq) for seq in sequences)
    padded_sequences = [seq.ljust(max_length, '-') for seq in sequences]
    sequences_matrix = np.array([list(seq) for seq in padded_sequences])

    return sequence_ids, sequences_matrix


def print_alignment(sequence_ids, alignment, output_file="aligned_sequences.fasta"):
    """
    Writes the aligned sequences to a FASTA file, preserving the original sequence IDs.

    Parameters:
    - sequence_ids (list of str): The original sequence IDs.
    - alignment (np.ndarray): The array of aligned sequences.
    - output_file (str): Path to the output FASTA file.
    """

    seq_records = []
    for idx, seq in enumerate(alignment):
        # Convert array back to string and remove trailing gaps
        sequence_str = ''.join(seq)
        seq_record = SeqRecord(Seq(sequence_str), id=sequence_ids[idx], description="")
        seq_records.append(seq_record)

    SeqIO.write(seq_records, output_file, "fasta")


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
    parser.add_argument('--temperature', type=float, default=10.0, help='Initial temperature for simulated annealing')
    parser.add_argument('--cooling_rate', type=float, default=0.999, help='Cooling rate for simulated annealing')
    parser.add_argument('--min_temperature', type=float, default=0.01, help='Minimal temperature of the system')
    parser.add_argument('--max_no_changes', type=int, default=100, help='Max consecutive no change events before early stopping')
    parser.add_argument('--match_score', type=float, default=1.0, help='')
    parser.add_argument('--mismatch_score', type=float, default=0.0, help='')
    parser.add_argument('--gap_score', type=float, default=-1.0, help='')

    return parser.parse_args()


def plot_title(args):
    input_file_path = args.fasta_file
    quality_function = args.quality_function
    temp = args.temperature
    cooling_rate = args.cooling_rate
    min_temp = args.min_temperature

    # Extracting the filename from the input file path
    input_filename = os.path.basename(input_file_path)

    num_sequences = 0
    max_length = 0
    for record in SeqIO.parse(input_file_path, "fasta"):
        num_sequences += 1
        max_length = max(max_length, len(record.seq))

    # Constructing the title string
    title = (
        f"MSASA Input={input_filename}, Sequences={num_sequences}, Max Length={max_length} aa, "
        f"Quality Function={quality_function}, "
        f"Temp={temp}°C, Cooling Rate={cooling_rate}, Min Temp:={min_temp}°C"
    )
    return title


def plot_charts_from_log(log_file_path, plot_title, max_no_change_events):
    # Load the data
    data = pd.read_csv(log_file_path, sep="\t")
    data['Timestamp'] = pd.to_datetime(data['Timestamp'])
    data['Cumulative_Accepted'] = data['Accepted'].cumsum()
    data['Acceptance_Rate'] = data['Cumulative_Accepted'] / (data.index + 1)

    # Prepare the figure
    fig, axs = plt.subplots(nrows=4, ncols=3, figsize=(18, 18), sharex=False)
    plt.suptitle(plot_title, fontsize=16, va="bottom")

    # Cooling Schedule Visualization
    axs[0, 0].plot(data["Iteration"], data["Temperature"], label="Temperature")
    axs[0, 0].set_title('Temperature over Iteration')
    axs[0, 0].set_xlabel('Iteration')
    axs[0, 0].set_ylabel('Temperature')
    ax2 = axs[0, 0].twinx()

    ax2.plot(data["Iteration"], data["Total_Accepted"], label="Total Accepted", color='green')
    ax2.plot(data["Iteration"], data["Total_Rejected"], label="Total Rejected", color='red')
    ax2.set_ylabel('Counts', color='green')  # Set the label for the second y-axis
    lines, labels = axs[0, 0].get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc="best")
    axs[0, 0].tick_params(axis='y', colors='blue')  # Set y-axis tick color to blue
    ax2.tick_params(axis='y', colors='green')  # Set secondary y-axis tick color to green
    axs[0, 0].grid(True)

    # Score over Temperature
    axs[0, 1].plot(data['Temperature'], data['Current_Score'], label='Current Score')
    axs[0, 1].plot(data['Temperature'], data['Best_Score'], label='Best Score', color='green')
    axs[0, 1].invert_xaxis()  # Inverting the X-axis
    axs[0, 1].set_title('Scores over Temperature')
    axs[0, 1].set_xlabel('Temperature')
    axs[0, 1].set_ylabel('Scores')
    axs[0, 1].grid(True)
    axs[0, 1].legend(loc="best")

    # Score over Iteration
    axs[0, 2].plot(data["Iteration"], data["Current_Score"], label="Current Score")
    axs[0, 2].plot(data["Iteration"], data["Best_Score"], label="Best Score", color='green')
    axs[0, 2].set_title('Scores over Iteration')
    axs[0, 2].set_xlabel('Iteration')
    axs[0, 2].set_ylabel('Scores')
    axs[0, 2].grid(True)
    axs[0, 2].legend(loc="best")

    axs[1, 0].plot(data['Iteration'], data['Historical_Score_Improvement'], label='Historical Score Improvement', color='blue')
    axs[1, 0].set_title('Total improvement over Iteration')
    axs[1, 0].set_xlabel('Iteration')
    axs[1, 0].set_ylabel("Score diff")
    axs[1, 0].grid(True)

    axs[1, 1].plot(data['Temperature'], data['Historical_Score_Improvement'], label='Historical Score Improvement', color='blue')
    axs[1, 1].set_title("Total improvement over Temperature")
    axs[1, 1].invert_xaxis()  # Inverting the X-axis
    axs[1, 1].set_xlabel('Iteration')
    axs[1, 1].set_ylabel('Score diff')
    axs[1, 1].grid(True)

    # Score Over Time
    axs[1, 2].plot(data["Timestamp"], data["Historical_Score_Improvement"], label="Historical Score Improvement", color='blue')
    axs[1, 2].set_title("Total improvement over Timestamp")
    axs[1, 2].set_xlabel("Timestamp")
    axs[1, 2].set_ylabel("Score diff")
    axs[1, 2].xaxis.set_major_formatter(FuncFormatter(ms_formatter))
    axs[1, 2].grid(True)

    # Acceptance Rate Over Iterations
    axs[2, 0].plot(data['Iteration'], data['Acceptance_Rate'], label='Acceptance Rate', color='black')
    axs[2, 0].set_title('Acceptance Rate over Iteration')
    axs[2, 2].set_xlabel("Iteration")
    axs[2, 0].set_ylabel('Rate %')
    axs[2, 0].grid(True)

    # Histogram of Score Changes
    # Creating the histogram with different colors for different ranges
    data_score_change = data["Score_Change"]
    bins = np.histogram_bin_edges(data_score_change, bins=30)

    axs[2, 1].hist(data_score_change[data_score_change < 0], bins=bins, color="red", edgecolor="black")
    axs[2, 1].hist(data_score_change[data_score_change >= 0], bins=bins, color="green", edgecolor="black")
    axs[2, 1].set_title("Histogram of Score Changes")
    axs[2, 1].set_xlabel("Score Change")
    axs[2, 1].set_ylabel("Frequency")
    axs[2, 1].grid(True)

    # Iteration Over Timestamp
    axs[2, 2].plot(data["Timestamp"], data["Iteration"], label="Iteration")
    axs[2, 2].plot(data["Timestamp"], data["No_Changes"], label="No-Changes events")
    axs[2, 2].set_title('Iteration Over Timestamp')
    axs[2, 2].set_xlabel('Timestamp')
    axs[2, 2].set_ylabel('Iteration')
    axs[2, 2].axhline(y=max_no_change_events, color="red", linestyle="--", label=f"Threshold {max_no_change_events}")
    axs[2, 2].xaxis.set_major_formatter(FuncFormatter(ms_formatter))
    axs[2, 2].grid(True)
    axs[2, 2].legend(loc="best")
    fig.autofmt_xdate()

    #
    axs[3, 0].plot(data['Iteration'], data['Iteration_Score'], label='Iteration Score')
    axs[3, 0].plot(data["Iteration"], data["New_Score"], label="New Score")
    axs[3, 0].plot(data["Iteration"], data["Best_Score"], label="Best Score", color="green")
    axs[3, 0].set_title('Score over Iteration')
    axs[3, 0].set_xlabel('Iteration')
    axs[3, 0].set_ylabel('Scores')
    axs[3, 0].grid(True)
    axs[3, 0].legend(loc="best")

    # Iteration/New scores over Temperature
    axs[3, 1].plot(data['Temperature'], data['Iteration_Score'], label='Iteration Score')
    axs[3, 1].plot(data["Temperature"], data["New_Score"], label="New Score")
    axs[3, 1].plot(data["Temperature"], data["Best_Score"], label="Best Score", color='green')
    axs[3, 1].invert_xaxis()  # Inverting the X-axis
    axs[3, 1].set_title('Score over Temperature')
    axs[3, 1].set_xlabel("Temperature")
    axs[3, 1].set_ylabel('Scores')
    axs[3, 1].grid(True)
    axs[3, 1].legend(loc="best")

    axs[3, 2].plot(data["Timestamp"], data["Iteration_Score"], label="Iteration Score")
    axs[3, 2].plot(data["Timestamp"], data["New_Score"], label="New Score")
    axs[3, 2].plot(data["Timestamp"], data["Best_Score"], label="Best Score", color='green')
    axs[3, 2].set_title('Score over Timestamp')
    axs[3, 2].set_xlabel('Timestamp')
    axs[3, 2].set_ylabel('Scores')
    axs[3, 2].xaxis.set_major_formatter(FuncFormatter(ms_formatter))
    axs[3, 2].grid(True)
    axs[3, 2].legend(loc="best")
    fig.autofmt_xdate()

    for ax in axs.flat:
        ax.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    # Adjust layout and save
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4, wspace=0.4, top=0.88)

    plt.savefig(log_file_path + ".png", dpi=300)
    plt.close()


def log_header():
    # Get the root logger
    logger = logging.getLogger()

    # Temporarily change the formatter to exclude timestamp
    original_formatter = logger.handlers[0].formatter
    no_timestamp_formatter = logging.Formatter("%(message)s")
    logger.handlers[0].setFormatter(no_timestamp_formatter)

    # Log the header
    logging.info(
        "Timestamp\tIteration\tTemperature\tCurrent_Score\tIteration_Score\tNew_Score\tScore_Change\tBest_Score\tHistorical_Score_Improvement\tTotal_Accepted\tTotal_Rejected\tAccepted\tNo_Changes"
    )

    # Revert the formatter to include timestamps for subsequent log entries
    logger.handlers[0].setFormatter(original_formatter)


def parse_fasta(fasta_file):
    """Parses a FASTA file and returns a dictionary of sequences."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def assert_sequences(original_sequences, aligned_sequences):
    # Assert that all aligned sequences have the same length
    aligned_lengths = [len(seq) for seq in aligned_sequences.values()]
    assert (
        len(set(aligned_lengths)) == 1
    ), "Not all aligned sequences have the same length."

    # Assert that the length of each original sequence matches its aligned version, ignoring gaps
    for seq_id, original_seq in original_sequences.items():
        aligned_seq = aligned_sequences.get(seq_id)
        assert aligned_seq is not None, f"Aligned sequence for {seq_id} not found."

        original_length = len(original_seq)
        aligned_length = len(aligned_seq.replace("-", ""))
        assert (
            original_length == aligned_length
        ), f"Length mismatch for {seq_id}: original ({original_length}) vs aligned ({aligned_length})"


def main():
    args = parse_arguments()

    logging.basicConfig(
        filename=args.log_file,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s.%(msecs)03d\t%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    log_header()

    sequence_ids, sequences = load_fasta(args.fasta_file)
    sa = SimulatedAnnealing(
        sequences,
        sequence_ids,
        extend=args.extend,
        temp=args.temperature,
        cooling_rate=args.cooling_rate,
        quality_function=args.quality_function,
        min_temp=args.min_temperature,
        no_changes_limit=args.max_no_changes,
        match_score=args.match_score,
        mismatch_score=args.mismatch_score,
        gap_score=args.gap_score,
    )

    sa.anneal()

    # Save aligned output in FASTA format
    print_alignment(sa.sequence_ids, sa.sequences, args.aligned_fasta_file)

    original_sequences = parse_fasta(args.fasta_file)
    aligned_sequences = parse_fasta(args.aligned_fasta_file)

    assert_sequences(original_sequences, aligned_sequences)

    # Create plots
    plot_charts_from_log(args.log_file, plot_title(args), args.max_no_changes)

    # Read the alignment from the aligned FASTA file
    alignment = AlignIO.read(args.aligned_fasta_file, "fasta")

    # Write the alignment to an MSF file
    # AlignIO.write(alignment, args.aligned_fasta_file + ".msf", "msf")
    AlignIO.write(alignment, args.aligned_fasta_file + ".clustal", "clustal")


if __name__ == "__main__":
    main()
