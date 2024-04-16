from itertools import combinations_with_replacement
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Tuple
import numpy as np
import random
import math
import logging
from Bio.Align import substitution_matrices
from .objective_functions import (
    ObjectiveFunction,
    Identity,
    Coincidences,
    Similarity,
    GlobalLocalAlignmentQuality,
)

class SimulatedAnnealing:
    MATCH_SCORE = 1
    MISMATCH_SCORE = -1
    GAP_PENALTY = -1

    def __init__(
        self,
        sequences,
        logger,
        extend=False,
        temp=100.0,
        cooling_rate=0.95,
        min_temp=0.001,
        quality_function="identity",
        no_changes_limit = 100,
        match_score = 1.0,
        mismatch_score = 0.0,
        gap_penalty = -1.0,
        changes = 1,
        iteration_neighbors = 1,
    ):
        self.sequences = sequences
        self.extend_sequences = extend
        self.temp = temp
        self.cooling_rate = cooling_rate
        self.min_temp = min_temp
        self.no_changes_limit = no_changes_limit
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        self.quality_function = quality_function
        self.logger = logger
        self.gap_char = ObjectiveFunction.GAP_CHAR_NP_ARRAY
        self.changes = changes
        self.iteration_neighbors = iteration_neighbors
        self.padding = {
            0: lambda seq, pad_length: np.pad(
                seq,
                (0, pad_length - len(seq)),
                constant_values=ObjectiveFunction.GAP_CHAR_NP_ARRAY,
            ),
            1: lambda seq, pad_length: np.pad(
                seq,
                (pad_length - len(seq), 0),
                constant_values=ObjectiveFunction.GAP_CHAR_NP_ARRAY,
            ),
        }

    def build_matrix(self) -> Dict[Tuple[str, str], float]:
        matches_dict = {}
        for res1_char, res2_char in combinations_with_replacement("-ARNDCQEGHILKMFPSTWYVBZX*", r=2):
            res1 = res1_char.encode()
            res2 = res2_char.encode()

            if res1 == ObjectiveFunction.GAP_SYMBOL and res2 == ObjectiveFunction.GAP_SYMBOL:
                matches_dict[(res1, res2)] = self.gap_penalty
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]
            elif res1 == ObjectiveFunction.GAP_SYMBOL or res2 == ObjectiveFunction.GAP_SYMBOL:
                matches_dict[(res1, res2)] = self.gap_penalty
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]
            elif res1 == res2:
                matches_dict[(res1, res2)] = self.match_score
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]
            else:
                matches_dict[(res1, res2)] = self.mismatch_score
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]

        return matches_dict

    def build_substitution_matrix(self, score_matrix, multiplier=1.0) -> Dict[Tuple[str, str], float]:
        matches_dict = {}
        for res1_char, res2_char in combinations_with_replacement("-ARNDCQEGHILKMFPSTWYVBZX*", r=2):
            res1 = res1_char.encode()
            res2 = res2_char.encode()

            if res1 == ObjectiveFunction.GAP_SYMBOL and res2 == ObjectiveFunction.GAP_SYMBOL:
                matches_dict[(res1, res2)] = self.gap_penalty
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]
            elif res1 == ObjectiveFunction.GAP_SYMBOL or res2 == ObjectiveFunction.GAP_SYMBOL:
                matches_dict[(res1, res2)] = self.gap_penalty
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]
            else:
                matches_dict[(res1, res2)] = multiplier * score_matrix[(chr(ord(res1)), chr(ord(res2)))]
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]

        return matches_dict

    def create_quality_instance(self):
        if self.quality_function == "identity":
            return Identity(match_score=self.match_score, mismatch_score=self.mismatch_score, gap_penalty=self.gap_penalty)
        elif self.quality_function == "coincidences":
            return Coincidences(match_score=self.match_score, mismatch_score=self.mismatch_score, gap_penalty=self.gap_penalty)
        elif self.quality_function == "similarity_blosum62":
            blosum62 = substitution_matrices.load("BLOSUM62")
            return Similarity(self.build_substitution_matrix(blosum62))
        elif self.quality_function == "similarity_pam250":
            pam250 = substitution_matrices.load("PAM250")
            return Similarity(self.build_substitution_matrix(pam250, -1.0))
        elif self.quality_function == "global":
            matches_matrix = self.build_matrix()
            return GlobalLocalAlignmentQuality(matches_matrix)
        elif self.quality_function == "local":
            matches_matrix = self.build_matrix()
            return GlobalLocalAlignmentQuality(matches_matrix)
        else:
            raise ValueError("Unknown quality function")

    def log_async(
        self,
        iteration: int,
        temp: float,
        sequence_length: int,
        current_score: float,
        iteration_score: float,
        new_score: float,
        score_change: float,
        initial_score: float,
        best_score: float,
        total_accepted: int,
        total_rejected: int,
        accepted: bool,
        no_changes: int,
        sequence_energy_hits: int,
        column_energy_hits: int,
    ):
        try:
            historical_score_improvement: float = best_score - initial_score
            acceptance = np.exp(score_change / temp)

            if initial_score != 0.0:
                historical_improvement_percentage: float = (historical_score_improvement / abs(initial_score)) * 100
                current_improvement_percentage: float = (score_change / abs(initial_score)) * 100
            else:
                historical_improvement_percentage: float = 0
                current_improvement_percentage: float = 0

            row_values = f"{iteration}\t{sequence_length}\t{temp:.7f}\t{current_score:+.7f}\t{iteration_score:+.7f}\t{new_score:+.7f}\t{score_change:+.7f}\t{current_improvement_percentage:+.4f}\t{best_score:+.7f}\t{historical_score_improvement:+.7f}\t{historical_improvement_percentage:+.4f}\t{total_accepted}\t{total_rejected}\t{accepted}\t{no_changes}\t{acceptance}\t{sequence_energy_hits}\t{column_energy_hits}"
            self.logger.info(row_values)
        except:
            pass

    def anneal(self):
        self.quality_function_instance = self.create_quality_instance()

        current_temp, iteration, no_changes, total_accepted, total_rejected = self.temp, 0, 0, 0, 0

        current_score_columns = self.sequences.transpose()

        self.current_score: float = self.quality_function_instance.energy(current_score_columns.tobytes(), self.sequences.shape[0], self.sequences.shape[1])

        initial_score, best_score = self.current_score, self.current_score

        with ThreadPoolExecutor(max_workers=1) as executor:  # Using one worker for logging
            while (current_temp := current_temp * self.cooling_rate) > self.min_temp and no_changes < self.no_changes_limit:
                iteration_score: float = self.current_score

                # Find the best new sequences:
                neighbors_best_score = None
                for i in np.arange(self.iteration_neighbors):
                    new_sequences_i_j = self.sequences.view()

                    for j in np.arange(np.random.randint(1, self.changes)):
                        new_sequences_i_j = self.generate_new_sequences(new_sequences_i_j, addition_deletion_prob=0.4)

                    new_sequences_columns = new_sequences_i_j.transpose()
                    new_score_i_j: float = self.quality_function_instance.energy(new_sequences_columns.tobytes(), new_sequences_i_j.shape[0], new_sequences_i_j.shape[1])

                    if (neighbors_best_score is None) or (new_score_i_j > neighbors_best_score):
                        new_sequences = new_sequences_i_j
                        new_score = new_score_i_j

                # Once the best new state is found, we proceed with the regular SA flow

                score_change = new_score - iteration_score

                if (accepted := self.should_accept(score_change, current_temp)):
                    self.sequences     = new_sequences
                    self.current_score = new_score
                    total_accepted     += 1

                    if new_score > best_score:
                        best_score = new_score
                else:
                    total_rejected += 1

                if iteration_score == self.current_score:
                    no_changes += 1
                else:
                    no_changes = 0

                executor.submit(
                    self.log_async,
                    iteration,
                    current_temp,
                    self.sequences.shape[1],
                    self.current_score,
                    iteration_score,
                    new_score,
                    score_change,
                    initial_score,
                    best_score,
                    total_accepted,
                    total_rejected,
                    accepted,
                    no_changes,
                    self.quality_function_instance.sequence_energy_hits,
                    self.quality_function_instance.column_energy_hits,
                )

                iteration += 1

    def should_accept(self, delta, current_temp) -> bool:
        return (delta > 0) or (np.random.rand() < np.exp(delta / current_temp))

    def generate_new_sequences(self, sequences, addition_deletion_prob=0.5):
        direction = np.random.randint(0, 2)
        gap_positions = np.argwhere(sequences == self.gap_char)

        if np.random.rand() < addition_deletion_prob or len(gap_positions) == 0:
            seq_index = np.random.randint(len(sequences))
            col_index = np.random.randint(0, sequences[seq_index].size - 1)

            new_sequences_list = [
                np.insert(seq, col_index, self.gap_char)
                if current_seq_index == seq_index
                else self.padding[direction](seq.copy(), len(seq) + 1)
                for current_seq_index, seq in enumerate(sequences)
            ]
        else:
            seq_index, col_index = random.choice(gap_positions)

            new_sequences_list = [
                self.padding[direction](np.delete(seq, col_index), len(seq))
                if current_seq_index == seq_index else seq.copy()
                for current_seq_index, seq in enumerate(sequences)
            ]

        new_sequences = np.array(new_sequences_list, dtype="|S1")

        while np.all(new_sequences[:, 0] == self.gap_char) or np.all(new_sequences[:, -1] == self.gap_char):
            if np.all(new_sequences[:, 0] == self.gap_char):
                new_sequences = new_sequences[:, 1:]
            if np.all(new_sequences[:, -1] == self.gap_char):
                new_sequences = new_sequences[:, :-1]

        return new_sequences
