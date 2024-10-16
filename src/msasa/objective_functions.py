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
    GAP_SYMBOL = ord(b'-')
    GAP_CHAR = '-'
    GAP_CHAR_NP_ARRAY = np.array([GAP_SYMBOL], dtype=np.uint8)

    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    def compute_column(self, column) -> float:
        pass

    @lru_cache(maxsize=1024)
    def get_column_energy(self, column_hash) -> float:
        column = np.frombuffer(column_hash, dtype=np.uint8)

        return self.compute_column(column)

    def divider(self, rows, columns) -> float:
        return 1.0

    @lru_cache(maxsize=1024)
    def energy(self, columns_bytes, sequences, residues) -> float:
        results = np.sum([self.get_column_energy(columns_bytes[sequences*i:(sequences*i)+sequences]) for i in range(residues)])

        return results / self.divider(sequences, residues)


def coincidences_calc(char_counts, match_score, mismatch_score, gap_penalty, gap_symbol):
    highest_count = max([count for ch, count in char_counts])
    if highest_count == 1:
        return np.array([mismatch_score] * len(char_counts))

    return np.array(
        [
            (
                highest_count * gap_penalty
                if max_group_char == gap_symbol
                else highest_count * match_score
            )
            for max_group_char, max_group_repetitions in char_counts
            if max_group_repetitions == highest_count
        ]
    )

class Coincidences(ObjectiveFunction):

    def __init__(self, match_score: float = 1.0, mismatch_score:float = -0.0, gap_penalty: float = -1.0):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        super().__init__()

    def divider(self, rows, columns) -> float:
        return columns

    def compute_column(self, column) -> float:
        scores = coincidences_calc(
            list(Counter(column).items()),
            self.match_score,
            self.mismatch_score,
            self.gap_penalty,
            ObjectiveFunction.GAP_SYMBOL,
        )

        return scores.sum() / len(column)

    def __str__(self) -> str:
        return "Coincidences"


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

        first_element = column[0]
        if np.all(column == first_element):
            return self.match_score

        return self.mismatch_score

    def __str__(self) -> str:
        return "Identity"


class Similarity(ObjectiveFunction):
    def __init__(self, similarity_matrix) -> None:
        self.similarity_matrix = similarity_matrix
        super().__init__()

    def compute_column(self, column) -> float:
        pair_counts = Counter(combinations(column, 2))

        score = np.sum(self.similarity_matrix[pair] * count for pair, count in pair_counts.items())

        return score

    def __str__(self) -> str:
        return "Similarity"


class GlobalLocalAlignmentQuality(ObjectiveFunction):
    def __init__(self, score_matrix) -> None:
        self.matrix = score_matrix
        super().__init__()

    def divider(self, rows, columns) -> float:
        return columns

    def compute_column(self, column) -> float:
        pair_counts = Counter(combinations(column, 2))
        score = np.sum(self.matrix[pair] * count for pair, count in pair_counts.items())

        return score

    def __str__(self) -> str:
        return "GlobalLocal"
