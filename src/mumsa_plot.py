import sys
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

if __name__ == "__main__":
    input_file = sys.argv[1]
    print(input_file)
    output_file = sys.argv[2]
    print(output_file)

    with open(input_file, "r") as file:
        lines = file.readlines()[1:]  # Skip the header line

    # Parsing the data into lists
    alignments = []
    scores = []

    for line in lines:
        parts = line.rsplit(maxsplit=1)
        alignments.append(
            "B"
            + parts[0]
            .strip()
            .replace(".fasta", "")
            .replace(".aln", "")
            .replace(".tfa", "")
            .replace(".alignment", "")
        )
        scores.append(float(parts[1]))

    # Plotting the data
    plt.figure(figsize=(10, 12))
    bars = plt.barh(alignments, scores, color='skyblue')

    plt.xlabel('Overlap Score')
    plt.title('Overlap Scores to the reference alignment')
    plt.gca().invert_yaxis()
    plt.grid(axis='x', linestyle='--', alpha=0.6)

    for bar, alignment in zip(bars, alignments):
        if '.msasa.' in alignment:
            # Highlight label background if it contains '.msasa.'
            plt.gca().text(bar.get_width(), bar.get_y() + bar.get_height()/2, f'{bar.get_width():.6f}',
                        va='center', ha='left', color='black', fontsize=8, backgroundcolor='yellow')
        elif '.reference.' in alignment:
            plt.gca().text(bar.get_width(), bar.get_y() + bar.get_height()/2, f'{bar.get_width():.1f}',
                        va='center', ha='left', color='black', fontsize=8)
        else:
            # Regular label without highlighting
            plt.gca().text(bar.get_width(), bar.get_y() + bar.get_height()/2, f'{bar.get_width():.6f}',
                        va='center', ha='left', color='black', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
