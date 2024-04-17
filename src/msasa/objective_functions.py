from collections import Counter
from typing import Dict, Tuple
import numpy as np
from itertools import combinations
from methodtools import lru_cache
from abc import ABC, abstractmethod

# https://github.com/Cantalapiedra/msa_conservation_index?tab=readme-ov-file
# https://www.biostars.org/p/5067/#5076
# https://www.majordifferences.com/2014/02/difference-between-pam-and-blosum-matrix_1.html


class ObjectiveFunction(ABC):
    GAP_SYMBOL = b'-'
    GAP_CHAR = '-'
    GAP_CHAR_NP_ARRAY = np.array(["-"], dtype="|S1")

    def __init__(self) -> None:
        super().__init__()
        self.sequence_energy_hits = 0
        self.column_energy_hits = 0

    @abstractmethod
    def compute_column(self, column) -> float :
        pass

    @lru_cache(maxsize=819200)
    def get_column_energy(self, column_hash) -> float:
        self.column_energy_hits += 1
        # print("\tcolumn_energy", column_hash)

        column = np.frombuffer(column_hash, dtype="|S1")

        return self.compute_column(column)

    def divider(self, rows, columns) -> float:
        return 1.0

    @lru_cache(maxsize=409600)
    def energy(self, columns_bytes, sequences, residues) -> float:
        self.sequence_energy_hits += 1
        # print("energy", columns_bytes)

        results = np.sum([self.get_column_energy(columns_bytes[sequences*i:(sequences*i)+sequences]) for i in range(residues)])

        return results / self.divider(sequences, residues)


class Coincidences(ObjectiveFunction):

    def __init__(self, match_score: float = 1.0, mismatch_score:float = -0.0, gap_penalty: float = -1.0):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        super().__init__()

    def divider(self, rows, columns) -> float:
        return columns

    def compute_column(self, column) -> float:
        char_counts = Counter(column)

        highest_count = max(char_counts.values())
        if highest_count == 1:
            return self.mismatch_score

        scores = np.array(
        [
            self.gap_penalty if max_group_char == ObjectiveFunction.GAP_SYMBOL else self.match_score
            for max_group_char, max_group_repetitions in char_counts.items()
            if max_group_repetitions == highest_count
        ])

        return (scores * highest_count).sum() / len(column)


class Identity(ObjectiveFunction):
    def __init__(self, match_score=1, mismatch_score=0, gap_penalty=-1):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        super().__init__()

    def divider(self, rows, columns) -> float:
        return columns

    def compute_column(self, column) -> float:
        if np.all(column == ObjectiveFunction.GAP_SYMBOL):
            return self.gap_penalty
        elif np.all(column == column[0]):
            return self.match_score
        else:
            return self.mismatch_score


class Similarity(ObjectiveFunction):
    def __init__(self, similarity_matrix) -> None:
        self.similarity_matrix = similarity_matrix
        super().__init__()

    def compute_column(self, column) -> float:
        return np.array([self.similarity_matrix[(res1, res2)] for res1, res2 in combinations(column, 2)]).sum()


class GlobalLocalAlignmentQuality(ObjectiveFunction):
    def __init__(self, score_matrix) -> None:
        self.matrix = score_matrix
        super().__init__()

    def divider(self, rows, columns) -> float:
        return columns

    def compute_column(self, column) -> float:
        return np.array([self.matrix[(res1, res2)] for res1, res2 in combinations(column, 2)]).sum()
