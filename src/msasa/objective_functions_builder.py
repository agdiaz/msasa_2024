from itertools import combinations_with_replacement
from typing import Dict, Tuple
from Bio.Align import substitution_matrices
from .objective_functions import (
    ObjectiveFunction,
    Identity,
    Coincidences,
    Similarity,
    GlobalLocalAlignmentQuality,
)

class Builder:

    def __init__(self, match_score: float, mismatch_score: float, gap_penalty: float) -> None:
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        self.max_value = max(self.match_score, self.gap_penalty, self.mismatch_score)
        self.min_value = min(self.match_score, self.gap_penalty, self.mismatch_score)

    def build_matrix(self) -> Dict[Tuple[str, str], float]:
        matches_dict = {}
        scale = lambda n: (n - self.min_value) / (self.max_value - self.min_value)

        for res1_char, res2_char in combinations_with_replacement("-ARNDCQEGHILKMFPSTWYVBZX*", r=2):
            res1 = ord(res1_char)
            res2 = ord(res2_char)

            if res1 == ObjectiveFunction.GAP_SYMBOL and res2 == ObjectiveFunction.GAP_SYMBOL:
                matches_dict[(res1, res2)] = scale(self.gap_penalty)
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]
            elif res1 == ObjectiveFunction.GAP_SYMBOL or res2 == ObjectiveFunction.GAP_SYMBOL:
                matches_dict[(res1, res2)] = scale(self.gap_penalty)
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]
            elif res1 == res2:
                matches_dict[(res1, res2)] = scale(self.match_score)
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]
            else:
                matches_dict[(res1, res2)] = scale(self.mismatch_score)
                matches_dict[(res2, res1)] = matches_dict[(res1, res2)]

        return matches_dict

    def build_substitution_matrix(self, score_matrix) -> Dict[Tuple[str, str], float]:
        max_value = max(self.max_value, score_matrix.max())
        min_value = min(self.min_value, score_matrix.min())
        scale = lambda n: (n - min_value) / (max_value - min_value)

        matches_dict = {}
        gap_substitution_value = scale(self.gap_penalty)

        for res1_char, res2_char in combinations_with_replacement("-ARNDCQEGHILKMFPSTWYVBZX*", r=2):
            res1 = ord(res1_char)
            res2 = ord(res2_char)

            if res1 == ObjectiveFunction.GAP_SYMBOL and res2 == ObjectiveFunction.GAP_SYMBOL:
                matches_dict[(res1, res2)] = gap_substitution_value
                matches_dict[(res2, res1)] = gap_substitution_value
            elif res1 == ObjectiveFunction.GAP_SYMBOL or res2 == ObjectiveFunction.GAP_SYMBOL:
                matches_dict[(res1, res2)] = gap_substitution_value
                matches_dict[(res2, res1)] = gap_substitution_value
            else:
                substitution_value = scale(score_matrix.get((res1_char, res2_char), scale(self.mismatch_score)))

                matches_dict[(res1, res2)] = substitution_value
                matches_dict[(res2, res1)] = substitution_value

        return matches_dict

    def create_quality_instance(self, quality_function):
        if quality_function == "identity":
            return Identity(match_score=self.match_score, mismatch_score=self.mismatch_score, gap_penalty=self.gap_penalty)
        elif quality_function == "coincidences":
            return Coincidences(match_score=self.match_score, mismatch_score=self.mismatch_score, gap_penalty=self.gap_penalty)
        elif quality_function == "similarity_blosum62":
            blosum62 = substitution_matrices.load("BLOSUM62")
            return Similarity(self.build_substitution_matrix(blosum62))
        elif quality_function == "similarity_gonnet":
            gonnet92 = substitution_matrices.load("GONNET1992")
            return Similarity(self.build_substitution_matrix(gonnet92))
        elif quality_function == "similarity_pam250":
            pam250 = substitution_matrices.load("PAM250")
            return Similarity(self.build_substitution_matrix(pam250))
        elif quality_function == "global":
            matches_matrix = self.build_matrix()
            return GlobalLocalAlignmentQuality(matches_matrix)
        elif quality_function == "local":
            matches_matrix = self.build_matrix()
            return GlobalLocalAlignmentQuality(matches_matrix)
        else:
            raise ValueError("Unknown quality function")
