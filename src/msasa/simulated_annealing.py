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

    def build_matrix(self) -> Dict[Tuple[str, str], float]:
        matches_dict = {}
        for res1_char, res2_char in combinations_with_replacement("-ABCDEFGHIJKLMNOPQRSTUVWXYZ", r=2):
            res1 = res1_char.encode()
            res2 = res2_char.encode()

            if res1 == ObjectiveFunction.GAP_SYMBOL and res2 == ObjectiveFunction.GAP_SYMBOL:
                matches_dict[(res1, res2)] = self.gap_penalty * 2
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

    def create_quality_instance(self):
        if self.quality_function == "identity":
            return Identity(match_score=self.match_score, mismatch_score=self.mismatch_score, gap_penalty=self.gap_penalty)
        elif self.quality_function == "coincidences":
            return Coincidences(match_score=self.match_score, mismatch_score=self.mismatch_score, gap_penalty=self.gap_penalty)
        elif self.quality_function == "similarity_blosum62":
            blosum62 = substitution_matrices.load("BLOSUM62")
            return Similarity(blosum62, gap_penalty=self.gap_penalty)
        elif self.quality_function == "similarity_pam250":
            pam250 = substitution_matrices.load("PAM250")
            return Similarity(pam250, gap_penalty=self.gap_penalty, multiplier=-1.0)
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
        temp:float,
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
    ):
        historical_score_improvement: float = best_score - initial_score
        acceptance = math.exp(current_score / temp)

        if initial_score != 0.0:
            historical_improvement_percentage: float = (historical_score_improvement / abs(initial_score)) * 100
            current_improvement_percentage: float = (score_change / abs(initial_score)) * 100
        else:
            historical_improvement_percentage: float = 0
            current_improvement_percentage: float = 0

        row_values = f"{iteration}\t{sequence_length}\t{temp:.10f}\t{current_score:.10f}\t{iteration_score:.10f}\t{new_score:.10f}\t{score_change:.10f}\t{current_improvement_percentage:.6f}\t{best_score:.10f}\t{historical_score_improvement:.10f}\t{historical_improvement_percentage:.6f}\t{total_accepted}\t{total_rejected}\t{accepted}\t{no_changes}\t{acceptance:.10f}"
        self.logger.info(row_values)

    def anneal(self):
#         print("""
# MSASA - Simulated Annealing for MSA (Adrián Díaz, Gabriela Minetti)
# Initial temperature = {temp}; Cooling rate = {cooling_rate}
# Match score = {match}; Mismatch score = {mismatch}; GAP score = {gap}
# Consecutive No-Change limit = {no_change};
# Energy function = {energy_func}
#               """.format(temp=self.temp, cooling_rate=self.cooling_rate, match=self.match_score, mismatch=self.mismatch_score, gap=self.gap_penalty, no_change=self.no_changes_limit, energy_func=self.quality_function))

        self.quality_function_instance = self.create_quality_instance()

        iteration, no_changes, total_accepted, total_rejected = 0, 0, 0, 0

        self.current_score: float = self.quality_function_instance.energy(self.sequences)
        initial_score, best_score = self.current_score, self.current_score

        with ThreadPoolExecutor(max_workers=1) as executor:  # Using one worker for logging
            while self.temp > self.min_temp:
                iteration_score: float = self.current_score

                new_state = self.generate_new_state_nf()
                new_score: float = self.quality_function_instance.energy(new_state)

                score_change: float = new_score - self.current_score

                accepted = self.should_accept(score_change)
                if accepted:
                    self.sequences = new_state.copy()
                    self.current_score = new_score
                    total_accepted += 1

                    if new_score > best_score:
                        best_score = new_score
                else:
                    total_rejected += 1

                if iteration_score == self.current_score:
                    no_changes += 1
                else:
                    no_changes = 0

                executor.submit(self.log_async, iteration, self.temp, self.sequences.shape[1], self.current_score, iteration_score, new_score, score_change, initial_score, best_score, total_accepted, total_rejected, accepted, no_changes,)

                self.temp *= self.cooling_rate
                iteration += 1

                if no_changes > self.no_changes_limit:
                    break

    def should_accept(self, delta) -> bool:
        return delta > 0 or random.random() < math.exp(delta / self.temp)


    def generate_new_state_nf(self):
        new_state = self.sequences

        gap_char = ObjectiveFunction.GAP_CHAR_NP_ARRAY
        gap_positions = np.argwhere(new_state == gap_char)

        action = 0 if random.random() > 0.5 else 1
        direction = 0 if random.random() > 0.5 else 1
        padding = {
            0: lambda seq, pad_length: np.pad(seq, (0, pad_length - len(seq)), constant_values=gap_char),
            1: lambda seq, pad_length: np.pad(seq, (pad_length - len(seq), 0), constant_values=gap_char),
        }

        if action == 0 or len(gap_positions) == 0:
            seq_index = random.randint(0, len(new_state) - 1)
            col_index = random.randint(0, new_state[seq_index].size - 1)

            new_sequences_list = []
            for seq_index_iterator in range(len(new_state)):
                if seq_index_iterator == seq_index:
                    new_sequences_list.append(np.insert(new_state[seq_index].copy(), col_index, gap_char))
                else:
                    new_sequences_list.append(padding[direction](new_state[seq_index_iterator].copy(), len(new_state[seq_index_iterator]) + 1))

            new_state = np.array(new_sequences_list)

        elif action == 1 and len(gap_positions) > 0:
            seq_index, col_index = random.choice([(s, c) for s, c in gap_positions if c < (new_state.shape[1] - 1)])

            new_sequences_list = []
            for seq_index_iterator in range(len(new_state)):
                if seq_index_iterator == seq_index:
                    new_sequences_list.append(padding[direction](np.delete(new_state[seq_index_iterator].copy(), col_index), len(new_state[seq_index])))
                else:
                    new_sequences_list.append(new_state[seq_index_iterator].copy())

            new_state = np.array(new_sequences_list)

        while np.all(new_state[:, 0] == gap_char) or np.all(new_state[:, -1] == gap_char):
            if np.all(new_state[:, 0] == gap_char):
                new_state = new_state[:, 1:]
            if np.all(new_state[:, -1] == gap_char):
                new_state = new_state[:, :-1]

        return new_state
