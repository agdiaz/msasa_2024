import os
import numpy as np
from Bio import SeqIO
from Bio import SeqIO

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt

import pandas as pd
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter
import seaborn as sns

def ms_formatter(a, b):
    t = mdates.num2date(a)
    ms = str(t.microsecond)[:3]
    res = f"{t.hour:02}:{t.minute:02}:{t.second:02}.{ms}"
    return res


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
    data["Timestamp"] = pd.to_datetime(data["Timestamp"])
    data["Cumulative_Accepted"] = data["Accepted"].cumsum()
    data["Acceptance_Rate"] = data["Cumulative_Accepted"] / (data.index + 1)

    # Prepare the figure
    fig, axs = plt.subplots(nrows=4, ncols=3, figsize=(18, 18), sharex=False)
    plt.suptitle(plot_title, fontsize=16, va="bottom")

    # Cooling Schedule Visualization
    axs[0, 0].plot(data["Iteration"], data["Temperature"], label="Temperature")
    axs[0, 0].set_title("Temperature over Iteration")
    axs[0, 0].set_xlabel("Iteration")
    axs[0, 0].set_ylabel("Temperature")
    ax2 = axs[0, 0].twinx()

    ax2.plot(
        data["Iteration"],
        data["Max_Length"],
        label="length",
        color="green",
    )
    ax2.set_ylabel("Length")  # Set the label for the second y-axis
    lines, labels = axs[0, 0].get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc="best")
    axs[0, 0].tick_params(axis="y")  # Set y-axis tick color to blue
    ax2.tick_params(axis="y")
    axs[0, 0].grid(True)

    # Score over Temperature
    axs[0, 1].plot(data["Temperature"], data["Current_Score"], label="Current Score")
    axs[0, 1].plot(
        data["Temperature"], data["Best_Score"], label="Best Score", color="green"
    )
    axs[0, 1].invert_xaxis()  # Inverting the X-axis
    axs[0, 1].set_title("Scores over Temperature")
    axs[0, 1].set_xlabel("Temperature")
    axs[0, 1].set_ylabel("Scores")
    axs[0, 1].grid(True)
    axs[0, 1].legend(loc="best")

    # Score over Iteration
    axs[0, 2].plot(data["Iteration"], data["Current_Score"], label="Current Score")
    axs[0, 2].plot(
        data["Iteration"], data["Best_Score"], label="Best Score", color="green"
    )
    axs[0, 2].set_title("Scores over Iteration")
    axs[0, 2].set_xlabel("Iteration")
    axs[0, 2].set_ylabel("Scores")
    axs[0, 2].grid(True)
    axs[0, 2].legend(loc="best")

    axs[1, 0].plot(
        data["Iteration"],
        data["Historical_Score_Improvement"],
        label="Historical Score Improvement",
        color="blue",
    )
    axs[1, 0].set_title("Total improvement over Iteration")
    axs[1, 0].set_xlabel("Iteration")
    axs[1, 0].set_ylabel("Score diff")
    axs[1, 0].grid(True)

    axs[1, 1].plot(
        data["Temperature"],
        data["Historical_Score_Improvement"],
        label="Historical Score Improvement",
        color="blue",
    )
    axs[1, 1].set_title("Total improvement over Temperature")
    axs[1, 1].invert_xaxis()  # Inverting the X-axis
    axs[1, 1].set_xlabel("Iteration")
    axs[1, 1].set_ylabel("Score diff")
    axs[1, 1].grid(True)

    #
    axs[1, 2].plot(data["Temperature"], data["Acceptance"], label="Acceptance Rate", color="black")
    axs[1, 2].set_title("Acceptance over Temperature")
    axs[1, 2].set_xlabel("Temperature")
    axs[1, 2].set_ylabel("Acceptance %")
    axs[1, 2].grid(True)

    # Acceptance Rate Over Iterations
    axs[2, 0].plot(data["Iteration"], data["Acceptance"], label="Acceptance Rate", color="black")
    axs[2, 0].set_title("Acceptance over Iteration")
    axs[2, 0].set_xlabel("Iteration")
    axs[2, 0].set_ylabel("Acceptance %")
    axs[2, 0].grid(True)

    # Histogram of Score Changes
    # Creating the histogram with different colors for different ranges
    data_score_change = data["Score_Change"]
    bins = np.histogram_bin_edges(data_score_change, bins=30)

    axs[2, 1].hist(
        data_score_change[data_score_change < 0],
        bins=bins,
        color="red",
        edgecolor="black",
    )
    axs[2, 1].hist(
        data_score_change[data_score_change >= 0],
        bins=bins,
        color="green",
        edgecolor="black",
    )
    axs[2, 1].set_title("Histogram of Score Changes")
    axs[2, 1].set_xlabel("Score Change")
    axs[2, 1].set_ylabel("Frequency")
    axs[2, 1].grid(True)

    # Iteration Over Timestamp
    axs[2, 2].plot(data["Timestamp"], data["Iteration"], label="Iteration")
    axs[2, 2].plot(data["Timestamp"], data["No_Changes"], label="No-Changes events")
    axs[2, 2].set_title("Iteration Over Timestamp")
    axs[2, 2].set_xlabel("Timestamp")
    axs[2, 2].set_ylabel("Iteration")
    axs[2, 2].axhline(
        y=max_no_change_events,
        color="red",
        linestyle="--",
        label=f"Threshold {max_no_change_events}",
    )
    axs[2, 2].xaxis.set_major_formatter(FuncFormatter(ms_formatter))
    axs[2, 2].grid(True)
    axs[2, 2].legend(loc="best")
    fig.autofmt_xdate()

    #
    axs[3, 0].plot(data["Iteration"], data["Iteration_Score"], label="Iteration Score")
    axs[3, 0].plot(data["Iteration"], data["New_Score"], label="New Score")
    axs[3, 0].plot(
        data["Iteration"], data["Best_Score"], label="Best Score", color="green"
    )
    axs[3, 0].set_title("Score over Iteration")
    axs[3, 0].set_xlabel("Iteration")
    axs[3, 0].set_ylabel("Scores")
    axs[3, 0].grid(True)
    axs[3, 0].legend(loc="best")

    # Iteration/New scores over Temperature
    axs[3, 1].plot(
        data["Temperature"], data["Iteration_Score"], label="Iteration Score"
    )
    axs[3, 1].plot(data["Temperature"], data["New_Score"], label="New Score")
    axs[3, 1].plot(
        data["Temperature"], data["Best_Score"], label="Best Score", color="green"
    )
    axs[3, 1].invert_xaxis()  # Inverting the X-axis
    axs[3, 1].set_title("Score over Temperature")
    axs[3, 1].set_xlabel("Temperature")
    axs[3, 1].set_ylabel("Scores")
    axs[3, 1].grid(True)
    axs[3, 1].legend(loc="best")

    axs[3, 2].plot(data["Timestamp"], data["Iteration_Score"], label="Iteration Score")
    axs[3, 2].plot(data["Timestamp"], data["New_Score"], label="New Score")
    axs[3, 2].plot(
        data["Timestamp"], data["Best_Score"], label="Best Score", color="green"
    )
    axs[3, 2].set_title("Score over Timestamp")
    axs[3, 2].set_xlabel("Timestamp")
    axs[3, 2].set_ylabel("Scores")
    axs[3, 2].xaxis.set_major_formatter(FuncFormatter(ms_formatter))
    axs[3, 2].grid(True)
    axs[3, 2].legend(loc="best")
    fig.autofmt_xdate()

    for ax in axs.flat:
        ax.tick_params(axis="x", which="both", bottom=True, top=False, labelbottom=True)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    # Adjust layout and save
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4, wspace=0.4, top=0.88)

    plt.savefig(log_file_path + ".png", dpi=300)
    plt.savefig(log_file_path + ".pdf", format='pdf', bbox_inches='tight')
    plt.close()

def plot_seaborn_charts_from_log(log_file_path, plot_title):
    data = pd.read_csv(log_file_path, sep="\t")
    data['Timestamp'] = pd.to_datetime(data['Timestamp'])

    # Create a large figure to hold all subplots
    fig = plt.figure(figsize=(18, 30))  # Adjusted size for readability
    plt.suptitle(plot_title, fontsize=16, va="bottom")

    # Scatter plot for Score over Iterations
    ax1 = fig.add_subplot(521)  # 5 rows, 2 columns, 1st subplot
    ax1.scatter(data['Iteration'], data['Best_Score'], color='blue', label='Best Score')
    ax1.scatter(data['Iteration'], data['Current_Score'], color='red', alpha=0.5, label='Current Score')
    ax1.set_title('Score Progression Over Iterations')
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Score')
    ax1.legend()
    ax1.grid(True)

    # Heat map of Temperature over Iterations
    ax2 = fig.add_subplot(522)  # 5 rows, 2 columns, 2nd subplot
    temperature_matrix = data.pivot(index='Iteration', columns='Timestamp', values='Temperature')
    sns.heatmap(temperature_matrix, cmap='coolwarm', ax=ax2)
    ax2.set_title('Temperature Change Over Iterations')
    ax2.set_xlabel('Timestamp')
    ax2.set_ylabel('Iteration')
    ax2.xaxis.set_major_formatter(FuncFormatter(ms_formatter))

    # Bar plot for Acceptance Rate per Iteration
    ax3 = fig.add_subplot(523)  # 5 rows, 2 columns, 3rd subplot
    ax3.bar(data['Iteration'], data['Acceptance'], color='green')
    ax3.set_title('Acceptance Rate per Iteration')
    ax3.set_xlabel('Iteration')
    ax3.set_ylabel('Acceptance Rate')

    # Line plot for Best Score over Time
    ax4 = fig.add_subplot(524)  # 5 rows, 2 columns, 4th subplot
    ax4.plot(data['Timestamp'], data['Best_Score'], marker='o', linestyle='-', color='magenta')
    ax4.set_title('Best Score Over Time')
    ax4.set_xlabel('Timestamp')
    ax4.xaxis.set_major_formatter(FuncFormatter(ms_formatter))
    ax4.set_xticklabels(ax4.get_xticklabels(), rotation=45)
    ax4.set_ylabel('Best Score')
    ax4.grid(True)

    # Temperature Decay
    ax5 = fig.add_subplot(525)  # 5 rows, 2 columns, 5th subplot
    ax5.plot(data['Iteration'], data['Temperature'], color='cyan')
    ax5.set_title('Temperature Decay Over Iterations')
    ax5.set_xlabel('Iteration')
    ax5.set_ylabel('Temperature')
    ax5.grid(True)

    # Histogram of Score Changes
    ax6 = fig.add_subplot(526)  # 5 rows, 2 columns, 6th subplot
    ax6.hist(data['Score_Change'], bins=30, color='purple', alpha=0.7)
    ax6.set_title('Histogram of Score Changes')
    ax6.set_xlabel('Score Change')
    ax6.set_ylabel('Frequency')
    ax6.grid(True)

    # Cumulative Count of Acceptances and Rejections
    ax7 = fig.add_subplot(527)  # 5 rows, 2 columns, 7th subplot
    ax7.plot(data['Iteration'], data['Total_Accepted'].cumsum(), label='Cumulative Acceptances', color='green')
    ax7.plot(data['Iteration'], data['Total_Rejected'].cumsum(), label='Cumulative Rejections', color='red')
    ax7.set_title('Cumulative Acceptances and Rejections Over Iterations')
    ax7.set_xlabel('Iteration')
    ax7.set_ylabel('Cumulative Count')
    ax7.legend()
    ax7.grid(True)

    # Box Plot for Scores by Iteration Blocks
    ax8 = fig.add_subplot(528)  # 5 rows, 2 columns, 8th subplot
    data['Iteration Block'] = pd.cut(data['Iteration'], bins=10)
    sns.boxplot(x='Iteration Block', y='Current_Score', data=data, ax=ax8)
    ax8.set_title('Box Plot of Scores by Iteration Blocks')
    ax8.set_xlabel('Iteration Block')
    ax8.set_ylabel('Current Score')
    ax8.set_xticklabels(ax8.get_xticklabels(), rotation=45)

    # Density Plot of Score Improvement Percentage
    ax9 = fig.add_subplot(529)  # 5 rows, 2 columns, 9th subplot
    sns.kdeplot(data['Score_Improvement_Perc'], shade=True, color='orange', ax=ax9)
    ax9.set_title('Density Plot of Score Improvement Percentage')
    ax9.set_xlabel('Score Improvement Percentage')
    ax9.set_ylabel('Density')

    # Adjust layout
    plt.tight_layout(pad=3.0)  # Adjust spacing to prevent overlap

    # Save the figure
    plt.savefig(log_file_path + ".sns.png", dpi=300)
    plt.savefig(log_file_path + ".sns.pdf", format='pdf', bbox_inches='tight')
    plt.close()