import numpy as np
import random
import math
import logging
from Bio.Align import substitution_matrices
from .objective_functions import (
    partial_identity,
    partial_coincidences,
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
        no_changes_limit = 100,
        match_score = 1.0,
        mismatch_score = 0.0,
        gap_score = -1.0
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
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score

        self.quality_function = self.select_quality_function(quality_function)

    def select_quality_function(self, quality_function_name):
        if quality_function_name == "identity":
            return partial_identity(match_score=self.match_score, mismatch_score=self.mismatch_score, gap_penalty=self.gap_score)
        elif quality_function_name == "coincidences":
            return partial_coincidences(match_score=self.match_score, gap_penalty=self.gap_score)
        elif quality_function_name == "similarity_blosum62":
            return partial_similarity(self.blosum62, gap_penalty=self.gap_score)
        elif quality_function_name == "similarity_pam250":
            return partial_similarity(self.pam250, gap_penalty=self.gap_score, multiplier=-1.0)
        elif quality_function_name == "global":
            return partial_global_alignment_quality(gap_penalty=self.gap_score, match_score=self.match_score, mismatch_score=self.mismatch_score)
        elif quality_function_name == "local":
            return partial_local_alignment_quality(gap_penalty=self.gap_score, match_score=self.match_score, mismatch_score=self.mismatch_score)
        else:
            raise ValueError("Unknown quality function")

    def anneal(self):
        iteration = 0
        no_changes = 0

        total_accepted = 0
        total_rejected = 0

        self.current_score = self.quality_function(self.sequences)
        initial_score = self.current_score
        best_score = self.current_score

        while self.temp > self.min_temp:
            iteration_score = self.current_score

            new_state = self.generate_new_state()
            new_score = self.quality_function(new_state)

            score_change = new_score - iteration_score

            if self.should_accept(new_score):
                self.sequences = new_state.copy()
                self.current_score = new_score

                accepted = True
                total_accepted += 1
                # score_changes.append(score_change)

                if new_score > best_score:
                    best_score = new_score
            else:
                accepted = False
                total_rejected += 1

            if iteration_score == self.current_score:
                no_changes += 1
            else:
                no_changes = 0

            # Log the required details for plotting
            historical_score_improvement = best_score - initial_score
            row_values = f"{iteration}\t{self.temp:.7f}\t{self.current_score:.7f}\t{iteration_score:.7f}\t{new_score:.7f}\t{score_change:.7f}\t{best_score:.7f}\t{historical_score_improvement:.7f}\t{total_accepted}\t{total_rejected}\t{accepted}\t{no_changes}"
            logging.info(row_values)

            self.temp *= self.cooling_rate
            iteration += 1

            if no_changes >= self.no_changes_limit:
                break

    def should_accept(self, new_score):
        if new_score >= self.current_score:
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
        all_same_length = all(len(seq) == max_length for seq in normalized_state)

        if not all_same_length:
            if random.random() > 0.5:
                normalized_state = [seq.ljust(max_length, '-') for seq in normalized_state]
            else:
                normalized_state = [seq.rjust(max_length, '-') for seq in normalized_state]

        should_trim = True
        while should_trim:
            # Normalize the sequences to ensure they all have the same length
            # max_length = max(len(seq) for seq in normalized_state)

            # Remove all-gaps-column at very beginning and very end of the sequences
            if all(seq[0] == '-' for seq in normalized_state):
                normalized_state = [seq[1:] for seq in normalized_state]
            elif all(seq[-1] == '-' for seq in normalized_state):
                normalized_state = [seq[:-1] for seq in normalized_state]
            else:
                should_trim = False

        return np.array([list(seq) for seq in normalized_state])
