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

    # Initialize DP table
    for i in range(n + 1):
        dp[i][0] = i * gap_penalty
    for j in range(m + 1):
        dp[0][j] = j * gap_penalty

    # Fill DP table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            char1, char2 = seq1[i - 1], seq2[j - 1]

            current_match_mismatch_score = 0
            if char1 == char2 and char1 != '-':  # Match
                current_match_mismatch_score = match_score
            elif char1 == '-' and char2 == '-':  # Both are gaps already (can happen if aligning profiles)
                current_match_mismatch_score = 0  # Or some other neutral/gap extension penalty
            elif char1 == '-' or char2 == '-':  # One is a gap, align with gap penalty
                current_match_mismatch_score = gap_penalty
            else:  # Mismatch
                current_match_mismatch_score = mismatch_penalty

            score_diag = dp[i - 1][j - 1] + current_match_mismatch_score
            score_up = dp[i - 1][j] + gap_penalty
            score_left = dp[i][j - 1] + gap_penalty
            dp[i][j] = max(score_diag, score_up, score_left)

    # Traceback
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
            else:  # Mismatch
                diag_path_score_component = mismatch_penalty

        # Check diagonal path (match/mismatch)
        if i_tb > 0 and j_tb > 0 and abs(
                current_score_val - (dp[i_tb - 1][j_tb - 1] + diag_path_score_component)) < 1e-9:
            aligned_s1.append(s1_char_tb)
            aligned_s2.append(s2_char_tb)
            i_tb -= 1
            j_tb -= 1
        # Check up path (gap in seq2)
        elif i_tb > 0 and abs(current_score_val - (dp[i_tb - 1][j_tb] + gap_penalty)) < 1e-9:
            aligned_s1.append(s1_char_tb)
            aligned_s2.append('-')
            i_tb -= 1
        # Check left path (gap in seq1)
        elif j_tb > 0 and abs(current_score_val - (dp[i_tb][j_tb - 1] + gap_penalty)) < 1e-9:
            aligned_s1.append('-')
            aligned_s2.append(s2_char_tb)
            j_tb -= 1
        else:  # Fallback should ideally not be common if logic is correct
            if i_tb > 0:  # Prefer consuming from seq1 if ambiguous
                aligned_s1.append(s1_char_tb)
                aligned_s2.append('-')
                i_tb -= 1
            elif j_tb > 0:
                aligned_s1.append('-')
                aligned_s2.append(s2_char_tb)
                j_tb -= 1
            else:  # Both i_tb and j_tb are 0
                break
    return score, "".join(reversed(aligned_s1)), "".join(reversed(aligned_s2))


def get_center_sequence_info(sequences_with_ids: list[tuple[str, str]], match_score: int = 1,
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
                pairwise_scores_matrix[i][j] = 0.0  # Score of a sequence with itself
                continue
            if j < i:  # Use already computed score
                score_ij = pairwise_scores_matrix[j][i]
            else:
                score_ij, _, _ = needleman_wunsch(seq_strings_only[i], seq_strings_only[j], match_score,
                                                  mismatch_penalty, gap_penalty)
            pairwise_scores_matrix[i][j] = score_ij
            if j > i: pairwise_scores_matrix[j][i] = score_ij  # Symmetric matrix
            current_i_total_score += score_ij
        total_scores_for_center_candidate[i] = current_i_total_score
        if current_i_total_score > max_total_score:
            max_total_score = current_i_total_score
            center_idx = i

    if center_idx == -1 and num_sequences > 0:  # Should always find a center if sequences exist
        center_idx = 0  # Default to first sequence if scores are all equal or problematic

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
        return [original_sequences_with_ids[0][1]]  # Return the single sequence as is

    # 1. Determine the center sequence
    center_original_idx_internal, center_item_internal, other_items_internal, _, _, _ = \
        get_center_sequence_info(original_sequences_with_ids, match_score, mismatch_penalty, gap_penalty)

    if center_item_internal is None:  # Should not happen if input is not empty
        return [""] * num_input_seqs

    # Initialize MSA with the center sequence (as a list of chars for easy modification)
    center_seq_str_current_msa = list(center_item_internal[1])
    msa_sequences_list = [center_seq_str_current_msa]  # List of lists of characters

    # 2. Align each other sequence to the current MSA (represented by the evolving center sequence/profile)
    for other_id, other_seq_raw in other_items_internal:
        center_profile_str = "".join(msa_sequences_list[0])  # Current state of the center sequence in MSA

        # Align the new sequence to the current center profile
        _, aligned_center_profile, aligned_new_other = needleman_wunsch(
            center_profile_str, other_seq_raw, match_score, mismatch_penalty, gap_penalty
        )

        # Update the MSA:
        # The center sequence in the MSA is updated first
        msa_sequences_list[0] = list(aligned_center_profile)

        # Propagate gaps from the new center profile to all other sequences already in MSA
        for k in range(1, len(msa_sequences_list)):  # For each existing non-center sequence in MSA
            current_other_seq_in_msa = msa_sequences_list[k]
            updated_other_seq = []
            old_profile_char_idx = 0  # Index for the center_profile_str (before this alignment)

            for new_profile_char in aligned_center_profile:  # Iterate through new center profile
                if old_profile_char_idx < len(center_profile_str) and \
                        new_profile_char == center_profile_str[old_profile_char_idx]:
                    # Character in center profile matches old profile char, copy from existing aligned sequence
                    updated_other_seq.append(current_other_seq_in_msa[old_profile_char_idx])
                    old_profile_char_idx += 1
                else:
                    # Gap introduced in center profile, or character changed (should be gap if new_profile_char != '-')
                    updated_other_seq.append('-')
                    # If the new_profile_char is NOT a gap, it means a character from other_seq_raw was aligned there.
                    # The old_profile_char_idx should still advance if the old char was not a gap.
                    if new_profile_char != '-' and \
                            old_profile_char_idx < len(center_profile_str) and \
                            center_profile_str[old_profile_char_idx] != '-':
                        old_profile_char_idx += 1

            msa_sequences_list[k] = updated_other_seq

        # Add the newly aligned sequence to the MSA
        msa_sequences_list.append(list(aligned_new_other))

    # 3. Reconstruct the final MSA strings in the original input order
    final_ordered_msa_strings = [""] * num_input_seqs
    if not msa_sequences_list or not msa_sequences_list[0]:  # MSA is empty
        return final_ordered_msa_strings

    # Place the (now aligned) center sequence in its original position
    final_ordered_msa_strings[center_original_idx_internal] = "".join(msa_sequences_list[0])

    # Place the other (now aligned) sequences in their original positions
    other_seq_msa_idx = 1  # Index for msa_sequences_list (starts from 1 as 0 is center)
    for i in range(num_input_seqs):
        if i == center_original_idx_internal:
            continue  # Skip center, already placed
        if other_seq_msa_idx < len(msa_sequences_list):
            final_ordered_msa_strings[i] = "".join(msa_sequences_list[other_seq_msa_idx])
            other_seq_msa_idx += 1
        else:
            # This case should ideally not be hit if logic is correct and all sequences processed
            # Fill with gaps if a sequence somehow got lost (highly unlikely)
            final_ordered_msa_strings[i] = "-" * len(msa_sequences_list[0]) if msa_sequences_list and \
                                                                               msa_sequences_list[0] else ""

    return final_ordered_msa_strings