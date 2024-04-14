from collections import Counter
from typing import Dict, Tuple
import numpy as np
from itertools import combinations
# from functools import lru_cache
from methodtools import lru_cache
from abc import ABC, abstractmethod

# https://github.com/Cantalapiedra/msa_conservation_index?tab=readme-ov-file
# https://www.biostars.org/p/5067/#5076
# https://www.majordifferences.com/2014/02/difference-between-pam-and-blosum-matrix_1.html


class ObjectiveFunction(ABC):
    GAP_SYMBOL = b'-'
    GAP_CHAR = '-'
    GAP_CHAR_NP_ARRAY = np.array(["-"], dtype="S1")

    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    @lru_cache(maxsize=128)
    def compute_column(self, column_hash) -> float:
        pass

    def divider(self, _sequences) -> float:
        return 1.0

    def energy(self, sequences) -> float:
        results = [self.compute_column(column.tobytes()) for column in sequences.transpose()]

        return sum(results) / self.divider(sequences)


class Coincidences(ObjectiveFunction):

    def __init__(self, match_score: float = 1.0, mismatch_score:float = -0.0, gap_penalty: float = -1.0):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty

    @lru_cache(maxsize=128)
    def compute_column(self, column_hash) -> float:
        column = np.frombuffer(column_hash, dtype="|S1")

        char_counts = Counter(column)

        highest_count = max(char_counts.values())
        if highest_count == 1:
            return self.mismatch_score

        score = 0

        for max_group_char, max_group_repetitions in char_counts.items():
            if max_group_repetitions == highest_count:
                if max_group_char == ObjectiveFunction.GAP_SYMBOL:
                    score += self.gap_penalty * highest_count
                else:
                    score += self.match_score * highest_count

        return score / len(column)


class Identity(ObjectiveFunction):
    def __init__(self, match_score=1, mismatch_score=0, gap_penalty=-1):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty

    def divider(self, sequences) -> float:
        return len(sequences)

    @lru_cache(maxsize=128)
    def compute_column(self, column_hash) -> float:
        column = np.frombuffer(column_hash, dtype="|S1")

        if np.all(column == ObjectiveFunction.GAP_SYMBOL):
            return self.gap_penalty
        elif np.all(column == column[0]):
            return self.match_score
        else:
            return self.mismatch_score


class Similarity(ObjectiveFunction):
    def __init__(self, similarity_matrix: Dict[Tuple[str, str], float], gap_penalty:float=-4.0, multiplier:float=1.0) -> None:
        self.similarity_matrix = similarity_matrix
        self.gap_penalty = gap_penalty
        self.multiplier = multiplier

    @lru_cache(maxsize=128)
    def score(self, res1: bytes, res2: bytes) -> float:
        if res1 == ObjectiveFunction.GAP_SYMBOL and res2 == ObjectiveFunction.GAP_SYMBOL:
            return self.gap_penalty * 2
        elif res1 == ObjectiveFunction.GAP_SYMBOL or res2 == ObjectiveFunction.GAP_SYMBOL:
            return self.gap_penalty
        else:
            return self.multiplier * self.similarity_matrix[(chr(ord(res1)), chr(ord(res2)))]

    @lru_cache(maxsize=128)
    def compute_column(self, column_hash) -> float:
        column = np.frombuffer(column_hash, dtype="|S1")

        return sum([self.score(res1, res2) for res1, res2 in combinations(column, 2)])


class GlobalLocalAlignmentQuality(ObjectiveFunction):
    def __init__(self, score_matrix: Dict[Tuple[str, str], float]) -> None:
        self.matrix = score_matrix

    @lru_cache(maxsize=128)
    def compute_column(self, column_hash) -> float:
        column = np.frombuffer(column_hash, dtype="|S1")
        return sum([self.matrix[(res1, res2)] for res1, res2 in combinations(column, 2)])
