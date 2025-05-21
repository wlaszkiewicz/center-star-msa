import numpy as np


def needleman_wunsch(seq1_str: str, seq2_str: str, match_score: int = 1, mismatch_penalty: int = -1,
                     gap_penalty: int = -1) -> tuple[float, str, str]:
    """
    Performs global alignment between two sequences using the Needleman-Wunsch algorithm.

    Args:
        seq1_str: The first sequence string.
        seq2_str: The second sequence string.
        match_score: Score for matching characters.
        mismatch_penalty: Penalty for mismatching characters.
        gap_penalty: Penalty for introducing a gap.

    Returns:
        A tuple containing:
            - The optimal alignment score (float).
            - The first aligned sequence string.
            - The second aligned sequence string.
    """
    seq1, seq2 = list(seq1_str), list(seq2_str)
    n, m = len(seq1), len(seq2)
    dp = np.zeros((n + 1, m + 1))

    for i in range(n + 1):
        dp[i][0] = i * gap_penalty
    for j in range(m + 1):
        dp[0][j] = j * gap_penalty

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            char1, char2 = seq1[i - 1], seq2[j - 1]

            current_match_mismatch_score = 0
            if char1 == char2 and char1 != '-':
                current_match_mismatch_score = match_score
            elif char1 == '-' and char2 == '-':
                current_match_mismatch_score = 0
            elif char1 == '-' or char2 == '-':
                current_match_mismatch_score = gap_penalty
            else:
                current_match_mismatch_score = mismatch_penalty

            score_diag = dp[i - 1][j - 1] + current_match_mismatch_score
            score_up = dp[i - 1][j] + gap_penalty
            score_left = dp[i][j - 1] + gap_penalty
            dp[i][j] = max(score_diag, score_up, score_left)

    aligned_s1, aligned_s2 = [], []
    i_tb, j_tb = n, m
    score = dp[n][m]

    while i_tb > 0 or j_tb > 0:
        current_score_val = dp[i_tb][j_tb]
        s1_char_tb = seq1[i_tb - 1] if i_tb > 0 else ''
        s2_char_tb = seq2[j_tb - 1] if j_tb > 0 else ''

        diag_path_score_component = 0
        if i_tb > 0 and j_tb > 0:
            if s1_char_tb == s2_char_tb and s1_char_tb != '-':
                diag_path_score_component = match_score
            elif s1_char_tb == '-' and s2_char_tb == '-':
                diag_path_score_component = 0
            elif s1_char_tb == '-' or s2_char_tb == '-':
                diag_path_score_component = gap_penalty
            else:
                diag_path_score_component = mismatch_penalty

        if i_tb > 0 and j_tb > 0 and abs(
                current_score_val - (dp[i_tb - 1][j_tb - 1] + diag_path_score_component)) < 1e-9:
            aligned_s1.append(s1_char_tb)
            aligned_s2.append(s2_char_tb)
            i_tb -= 1
            j_tb -= 1
        elif i_tb > 0 and abs(current_score_val - (dp[i_tb - 1][j_tb] + gap_penalty)) < 1e-9:
            aligned_s1.append(s1_char_tb)
            aligned_s2.append('-')
            i_tb -= 1
        elif j_tb > 0 and abs(current_score_val - (dp[i_tb][j_tb - 1] + gap_penalty)) < 1e-9:
            aligned_s1.append('-')
            aligned_s2.append(s2_char_tb)
            j_tb -= 1
        else:
            if i_tb > 0:
                aligned_s1.append(s1_char_tb)
                aligned_s2.append('-')
                i_tb -= 1
            elif j_tb > 0:
                aligned_s1.append('-')
                aligned_s2.append(s2_char_tb)
                j_tb -= 1
            else:
                break
    return score, "".join(reversed(aligned_s1)), "".join(reversed(aligned_s2))


def get_center_sequence(sequences_with_ids: list[tuple[str, str]], match_score: int = 1,
                        mismatch_penalty: int = -1, gap_penalty: int = -1) -> tuple[
    int, tuple[str, str] | None, list[tuple[str, str]], list[tuple[str, str]], list[list[float]], list[float]]:
    """
    Determines the center sequence from a list of sequences based on pairwise alignment scores.

    Args:
        sequences_with_ids: A list of tuples, where each tuple is (sequence_id, sequence_string).
        match_score: Score for matching characters in pairwise alignment.
        mismatch_penalty: Penalty for mismatching characters.
        gap_penalty: Penalty for introducing a gap.

    Returns:
        A tuple containing:
            - Index of the center sequence (int).
            - The center sequence tuple (id, sequence) or None if no sequences.
            - A list of other sequence tuples (excluding the center).
            - The original list of sequences_with_ids.
            - Pairwise alignment scores matrix (list of lists of float).
            - Total alignment scores for each sequence (list of float).
    """
    if not sequences_with_ids:
        return -1, None, [], sequences_with_ids, [[]], []

    num_sequences = len(sequences_with_ids)
    seq_strings_only = [item[1] for item in sequences_with_ids]

    if num_sequences == 1:
        return 0, sequences_with_ids[0], [], sequences_with_ids, [[0.0]], [0.0]

    pairwise_scores_matrix = [[0.0] * num_sequences for _ in range(num_sequences)]
    total_scores_for_center_candidate = [0.0] * num_sequences
    max_total_score = -float('inf')
    center_idx = -1

    for i in range(num_sequences):
        current_i_total_score = 0.0
        for j in range(num_sequences):
            if i == j:
                pairwise_scores_matrix[i][j] = 0.0
                continue
            if j < i:
                score_ij = pairwise_scores_matrix[j][i]
            else:
                score_ij, _, _ = needleman_wunsch(seq_strings_only[i], seq_strings_only[j], match_score,
                                                  mismatch_penalty, gap_penalty)
            pairwise_scores_matrix[i][j] = score_ij
            if j > i: pairwise_scores_matrix[j][i] = score_ij
            current_i_total_score += score_ij
        total_scores_for_center_candidate[i] = current_i_total_score
        if current_i_total_score > max_total_score:
            max_total_score = current_i_total_score
            center_idx = i

    if center_idx == -1 and num_sequences > 0:
        center_idx = 0

    center_item_tuple = sequences_with_ids[center_idx]
    other_items = [item for idx, item in enumerate(sequences_with_ids) if idx != center_idx]

    return center_idx, center_item_tuple, other_items, sequences_with_ids, pairwise_scores_matrix, total_scores_for_center_candidate


def center_star(original_sequences_with_ids: list[tuple[str, str]], match_score: int = 1, mismatch_penalty: int = -1,
                gap_penalty: int = -1) -> list[str]:
    """
    Performs multiple sequence alignment using the Center Star algorithm.

    Args:
        original_sequences_with_ids: List of (id, sequence) tuples.
        match_score: Score for matching characters.
        mismatch_penalty: Penalty for mismatching characters.
        gap_penalty: Penalty for introducing a gap.

    Returns:
        A list of aligned sequence strings, in the same order as input.
    """
    if not original_sequences_with_ids:
        return []

    num_input_seqs = len(original_sequences_with_ids)
    if num_input_seqs == 1:
        return [original_sequences_with_ids[0][1]]

    center_original_idx_internal, center_item_internal, other_items_internal, _, _, _ = \
        get_center_sequence(original_sequences_with_ids, match_score, mismatch_penalty, gap_penalty)

    if center_item_internal is None:
        return [""] * num_input_seqs

    center_seq_str_current_msa = list(center_item_internal[1])
    msa_sequences_list = [center_seq_str_current_msa]

    for other_id, other_seq_raw in other_items_internal:
        center_profile_str = "".join(msa_sequences_list[0])

        _, aligned_center_profile, aligned_new_other = needleman_wunsch(
            center_profile_str, other_seq_raw, match_score, mismatch_penalty, gap_penalty
        )

        msa_sequences_list[0] = list(aligned_center_profile)

        for k in range(1, len(msa_sequences_list)):
            current_other_seq_in_msa = msa_sequences_list[k]
            updated_other_seq = []
            old_profile_char_idx = 0

            for new_profile_char in aligned_center_profile:
                if old_profile_char_idx < len(center_profile_str) and \
                        new_profile_char == center_profile_str[old_profile_char_idx]:
                    updated_other_seq.append(current_other_seq_in_msa[old_profile_char_idx])
                    old_profile_char_idx += 1
                else:
                    updated_other_seq.append('-')
                    if new_profile_char != '-' and \
                            old_profile_char_idx < len(center_profile_str) and \
                            center_profile_str[old_profile_char_idx] != '-':
                        old_profile_char_idx += 1

            msa_sequences_list[k] = updated_other_seq

        msa_sequences_list.append(list(aligned_new_other))


    final_ordered_msa_strings = [""] * num_input_seqs
    if not msa_sequences_list or not msa_sequences_list[0]:
        return final_ordered_msa_strings

    final_ordered_msa_strings[center_original_idx_internal] = "".join(msa_sequences_list[0])

    other_seq_msa_idx = 1
    for i in range(num_input_seqs):
        if i == center_original_idx_internal:
            continue
        if other_seq_msa_idx < len(msa_sequences_list):
            final_ordered_msa_strings[i] = "".join(msa_sequences_list[other_seq_msa_idx])
            other_seq_msa_idx += 1
        else:
            final_ordered_msa_strings[i] = "-" * len(msa_sequences_list[0]) if msa_sequences_list and msa_sequences_list[0] else ""

    return final_ordered_msa_strings