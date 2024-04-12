from functools import partial
from collections import Counter
import numpy as np

# https://github.com/Cantalapiedra/msa_conservation_index?tab=readme-ov-file
# https://www.biostars.org/p/5067/#5076
# https://www.majordifferences.com/2014/02/difference-between-pam-and-blosum-matrix_1.html

def calculate_score(residue_a, residue_b, match_score, mismatch_score):
    if residue_a == residue_b:
        return match_score
    else:
        return mismatch_score


# def coincidences_per_column(state):
#     """
#     Counts the number of coincidences per column across all sequences.
#     """
#     coincidences = 0
#     for col in range(state.shape[1]):
#         column = state[:, col]
#         most_common = np.bincount(column).argmax()
#         coincidences += (column == most_common).sum()

#     return coincidences

# Checked !
def coincidences_per_column_with_ties_and_gap_penalty(state, gap_penalty=-1):
    """
    Counts the number of coincidences per column across all sequences in a numpy matrix of characters,
    including all characters that are tied for most common in their counts, and applies a gap penalty.

    Parameters:
    - state (np.ndarray): The state, a 2D numpy array of characters, including gaps.
    - gap_penalty (int): The penalty applied for each gap in a column.

    Returns:
    - int: The total coincidences score adjusted for gap penalties.
    """
    coincidences = 0
    for col in range(state.shape[1]):
        column = state[:, col]
        char_counts = Counter(column)

        # Find the highest frequency of characters in this column
        highest_count = max(char_counts.values())

        # Check if the most common character is a gap
        if char_counts.most_common(1)[0][0] == "-":
            # Apply the gap penalty for each gap in the column instead of adding coincidences
            coincidences += gap_penalty * highest_count
        else:
            # Count how many characters have this highest frequency, excluding gaps
            tied_chars_count = sum(count == highest_count for char, count in char_counts.items() if char != "-")

            # Add the highest count times the number of tied characters to the total coincidences
            coincidences += highest_count * tied_chars_count

    return coincidences

# Checked !
def identity(state, gap_penalty=-1):
    # Initialize counters
    identical_columns = 0
    total_columns = state.shape[1]

    # Iterate over each column
    for col in range(total_columns):
        column = state[:, col]
        # Check if all elements in the column are identical
        if np.all(column == "-"):
            identical_columns += gap_penalty
        elif np.all(column == column[0]):
            identical_columns += 1

    # Calculate identity percentage
    identity_percentage = (identical_columns / total_columns) * 100

    return identity_percentage


# def identity_with_gap_penalty(state, gap_penalty=-1):
#     # Initialize counters
#     identical_columns = 0
#     gap_columns = 0
#     total_columns = state.shape[1]

#     # Iterate over each column
#     for col in range(total_columns):
#         column = state[:, col]

#         # Check for columns with gaps
#         if "-" in column:
#             gap_columns += 1
#             # Additional condition to check if all are gaps, which could be treated differently
#             if np.all(column == "-"):
#                 continue  # or apply a different penalty

#         # Check if all non-gap elements in the column are identical
#         non_gap_elements = column[column != "-"]
#         if len(non_gap_elements) > 0 and np.all(non_gap_elements == non_gap_elements[0]):
#             identical_columns += 1

#     # Adjust identity calculation
#     adjusted_identical_columns = identical_columns + (gap_penalty * gap_columns)

#     # Calculate identity percentage
#     identity_percentage = (adjusted_identical_columns / total_columns) * 100

#     return identity_percentage


def partial_similarity(similarity_matrix, gap_penalty=-4):

    def similarity(state):
        similarity_score = 0

        for col in range(state.shape[1]):
            column = state[:, col]
            # Compare each pair of residues in the column
            for i in range(len(column)):
                for j in range(i + 1, len(column)):
                    if column[i] == "-" or column[j] == "-":
                        similarity_score += gap_penalty
                    else:
                        pair = (column[i], column[j])
                        similarity_score += similarity_matrix[
                            pair
                        ]  # Directly access the score

        return similarity_score

    return partial(similarity)


def partial_global_alignment_quality(gap_penalty, match_score, mismatch_score):

    def global_alignment_quality(state):
        quality_score = 0

        for col in range(state.shape[1]):
            column_scores = 0

            column = state[:, col]
            for i in range(len(column)):
                for j in range(i + 1, len(column)):
                    if "-" in (column[i], column[j]):
                        score = gap_penalty
                    else:
                        score = calculate_score(column[i], column[j], match_score, mismatch_score)

                    quality_score += score

            quality_score += column_scores

        return quality_score

    return partial(global_alignment_quality)


def partial_local_alignment_quality(match_score, mismatch_score, gap_penalty):

    def local_alignment_quality(state):
        quality_score = 0

        for col in range(state.shape[1]):
            column_scores = 0  # Collect positive scores for this column

            column = state[:, col]
            if np.all(column == "-"):
                quality_score += gap_penalty
                continue

            for i in range(len(column)):
                for j in range(i + 1, len(column)):
                    if "-" in (column[i], column[j]):
                        score = gap_penalty
                    else:
                        score = calculate_score(column[i], column[j], match_score, mismatch_score)

                    column_scores += score

            if column_scores > 0:
                quality_score += column_scores

        return quality_score

    return partial(local_alignment_quality)
