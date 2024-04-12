import numpy as np
import random
import math
import logging
from Bio.Align import substitution_matrices
from .objective_functions import (
    identity,
    coincidences_per_column_with_ties_and_gap_penalty,
    partial_similarity,
    partial_global_alignment_quality,
    partial_local_alignment_quality,
)

class SimulatedAnnealing:
    MATCH_SCORE = 1
    MISMATCH_SCORE = -1
    GAP_PENALTY = -1

    def __init__(
        self,
        sequences,
        sequence_ids,
        extend=False,
        temp=100.0,
        cooling_rate=0.95,
        min_temp=0.001,
        quality_function="identity",
        no_changes_limit = 100
    ):
        self.blosum62 = substitution_matrices.load("BLOSUM62")
        self.pam250 = substitution_matrices.load("PAM250")

        self.sequences = sequences
        self.sequence_ids = sequence_ids
        self.extend_sequences = extend
        self.temp = temp
        self.cooling_rate = cooling_rate
        self.min_temp = min_temp
        self.no_changes_limit = no_changes_limit
        self.quality_function = self.select_quality_function(quality_function)

    def select_quality_function(self, quality_function_name):
        if quality_function_name == "identity":
            return identity
        elif quality_function_name == "coincidences":
            return coincidences_per_column_with_ties_and_gap_penalty
        elif quality_function_name == "similarity_blosum62":
            return partial_similarity(self.blosum62, gap_penalty=0)
        elif quality_function_name == "similarity_pam250":
            return partial_similarity(self.pam250, gap_penalty=0)
        elif quality_function_name == "global":
            return partial_global_alignment_quality(gap_penalty=-1, match_score=1, mismatch_score=0)
        elif quality_function_name == "local":
            return partial_local_alignment_quality(gap_penalty=-1, match_score=2, mismatch_score=1)
        else:
            raise ValueError("Unknown quality function")

    def anneal(self):
        iteration = 0
        no_changes = 0

        total_accepted = 0
        # To keep track of score changes for histogram
        # score_changes = []

        self.current_score = self.quality_function(self.sequences)
        best_score = self.current_score

        while self.temp > self.min_temp:
            iteration_score = self.current_score

            new_state = self.generate_new_state()
            new_score = self.quality_function(new_state)

            score_change = new_score - iteration_score
            state_change_magnitude = self.calculate_state_change_magnitude(self.sequences, new_state)

            accepted = False
            if self.should_accept(new_score):
                self.sequences = new_state.copy()
                self.current_score = new_score

                accepted = True
                total_accepted += 1
                # score_changes.append(score_change)

                if new_score > best_score:
                    best_score = new_score

            if iteration_score == self.current_score:
                no_changes += 1
            else:
                no_changes = 0

            # Log the required details for plotting
            row_values = f"{iteration}\t{self.temp:.7f}\t{self.current_score:.7f}\t{iteration_score:.7f}\t{new_score:.7f}\t{best_score:.7f}\t{total_accepted}\t{score_change:.7f}\t{state_change_magnitude}\t{accepted}\t{no_changes}"
            logging.info(row_values)

            self.temp *= self.cooling_rate
            iteration += 1

            if no_changes >= self.no_changes_limit:
                break

    def calculate_state_change_magnitude(self, old_state, new_state):
        # Calculate overlaps in rows and columns to handle differing shapes
        row_overlap = min(old_state.shape[0], new_state.shape[0])
        col_overlap = min(old_state.shape[1], new_state.shape[1])

        # Count mismatches in the overlapping area
        mismatches = np.sum(old_state[:row_overlap, :col_overlap] != new_state[:row_overlap, :col_overlap])

        # For non-overlapping areas (if any), count all elements as mismatches
        # This includes extra rows or columns in either matrix
        if old_state.shape != new_state.shape:
            extra_rows_old = max(0, old_state.shape[0] - new_state.shape[0])
            extra_cols_old = max(0, old_state.shape[1] - new_state.shape[1])
            extra_rows_new = max(0, new_state.shape[0] - old_state.shape[0])
            extra_cols_new = max(0, new_state.shape[1] - old_state.shape[1])

            # Count extra elements as mismatches
            mismatches += extra_rows_old * old_state.shape[1] + extra_cols_old * old_state.shape[0]
            mismatches += extra_rows_new * new_state.shape[1] + extra_cols_new * new_state.shape[0]

            # Adjust for double-counting in the case where both matrices have extra rows/cols
            mismatches -= extra_rows_old * extra_cols_old
            mismatches -= extra_rows_new * extra_cols_new

        return mismatches

    def should_accept(self, new_score):
        if new_score > self.current_score:
            return True
        else:
            delta = new_score - self.current_score

            probability = math.exp(delta / self.temp)
            return random.random() < probability

    def generate_new_state(self):
        """
        Modifies the sequences by randomly adding or removing a gap in one sequence at a random position,
        then normalizes all sequences to the same length by adding gaps as needed.
        Assumes 'sequences' is a list of strings.
        """

        new_state = ["".join(seq) for seq in self.sequences]

        action = "add" if random.random() > 0.5 else "remove"
        # Collect positions of all gaps
        gap_positions = [(i, j) for i, seq in enumerate(new_state) for j, char in enumerate(seq) if char == "-"]

        if action == "add" or len(gap_positions) == 0:
            # Add a gap in a random position in one of the sequences
            seq_index = random.randint(0, len(new_state) - 1)
            col_index = random.randint(0, len(new_state[seq_index]))
            new_state[seq_index] = new_state[seq_index][:col_index] + "-" + new_state[seq_index][col_index:]
        elif action == "remove" and len(gap_positions) > 0:
            selected_gap = random.choice(gap_positions)
            seq_with_gap_removed = new_state[selected_gap[0]]
            new_state[selected_gap[0]] = seq_with_gap_removed[:selected_gap[1]] + seq_with_gap_removed[selected_gap[1] + 1:]

        normalized_state = new_state.copy()

        max_length = max(len(seq) for seq in normalized_state)
        if random.random() > 0.5:
            normalized_state = [seq.ljust(max_length, '-') for seq in normalized_state]
        else:
            normalized_state = [seq.rjust(max_length, '-') for seq in normalized_state]

        should_trim = True
        while should_trim:
            # Normalize the sequences to ensure they all have the same length
            max_length = max(len(seq) for seq in normalized_state)

            # Remove all-gaps-column at very beginning and very end of the sequences
            if all(seq[0] == '-' for seq in normalized_state):
                normalized_state = [seq[1:] for seq in normalized_state]
            elif all(seq[-1] == '-' for seq in normalized_state):
                normalized_state = [seq[:-1] for seq in normalized_state]
            else:
                should_trim = False

        # import ipdb; ipdb.set_trace()
        return np.array([list(seq) for seq in normalized_state])
