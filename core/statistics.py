def calculate_msa_statistics(aligned_sequences: list[str], match_score: int, mismatch_penalty: int, gap_penalty: int) -> \
tuple[float, int, int, float, list[list[float]], list[list[float]], int, int]:
    """
    Calculates various statistics for a given multiple sequence alignment.

    Args:
        aligned_sequences: A list of aligned sequence strings of equal length.
        match_score: Score used for a match in SP score calculation.
        mismatch_penalty: Penalty for a mismatch in SP score calculation.
        gap_penalty: Penalty for a gap in SP score calculation.

    Returns:
        A tuple containing:
            - Average pairwise identity (float, percentage).
            - Total number of gaps in the MSA (int).
            - Number of fully conserved columns (int).
            - Sum-of-Pairs (SP) score for the MSA (float).
            - Pairwise identity matrix (list of lists of float, percentage).
            - Pairwise distance matrix (list of lists of float, percentage).
            - Number of sequences in the alignment (int).
            - Length of the alignment (int).
    """
    if not aligned_sequences or not all(isinstance(s, str) for s in aligned_sequences):  # Basic check
        return 0.0, 0, 0, 0.0, [[]], [[]], 0, 0

    num_seqs = len(aligned_sequences)
    if num_seqs == 0:
        return 0.0, 0, 0, 0.0, [[]], [[]], 0, 0

    aln_len = len(aligned_sequences[0])
    if any(len(s) != aln_len for s in aligned_sequences):
        # This indicates a malformed MSA if sequences have different lengths
        # Return empty/error state or raise an exception
        return 0.0, 0, 0, 0.0, \
            [[0.0] * num_seqs for _ in range(num_seqs)], \
            [[0.0] * num_seqs for _ in range(num_seqs)], \
            num_seqs, 0  # Indicate 0 alignment length due to inconsistency

    if aln_len == 0:  # All sequences are empty strings
        empty_id_matrix = [[100.0 if i == j else 0.0 for j in range(num_seqs)] for i in range(num_seqs)]
        empty_dist_matrix = [[0.0 if i == j else 100.0 for j in range(num_seqs)] for i in range(num_seqs)]
        return 0.0, 0, 0, 0.0, empty_id_matrix, empty_dist_matrix, num_seqs, aln_len

    # 1. Total gaps in MSA
    total_gaps_msa = sum(s.count('-') for s in aligned_sequences)

    # 2. Fully conserved columns
    fully_conserved_cols = 0
    for j in range(aln_len):  # Iterate over columns
        first_char_in_col = None
        is_conserved_this_col = True
        has_any_nongap_char = False
        for i in range(num_seqs):  # Iterate over sequences for this column
            char_ij = aligned_sequences[i][j]
            if char_ij != '-':
                has_any_nongap_char = True
                if first_char_in_col is None:
                    first_char_in_col = char_ij
                elif first_char_in_col != char_ij:
                    is_conserved_this_col = False
                    break  # Column is not conserved
            # If char_ij is '-', it doesn't break conservation unless all are '-'
        if is_conserved_this_col and has_any_nongap_char:  # All non-gap chars are same
            fully_conserved_cols += 1

    # 3. Average Pairwise Identity (PID) & Matrices
    total_pid_sum = 0.0
    num_pairs = 0
    id_matrix = [[0.0] * num_seqs for _ in range(num_seqs)]
    dist_matrix = [[0.0] * num_seqs for _ in range(num_seqs)]

    for i in range(num_seqs):
        # Diagonal elements: identity is 100% (or 0% if seq is all gaps), distance is 0%
        if aln_len > 0 and aligned_sequences[i].count('-') < aln_len:  # Sequence is not all gaps
            id_matrix[i][i] = 100.0
        else:  # Sequence is all gaps or alignment length is 0
            id_matrix[i][i] = 0.0
        dist_matrix[i][i] = 0.0

        for j in range(i + 1, num_seqs):
            matches = 0
            relevant_cols_for_pid = 0  # Columns where at least one sequence is not a gap

            for k in range(aln_len):
                c1 = aligned_sequences[i][k]
                c2 = aligned_sequences[j][k]

                if c1 == '-' and c2 == '-':  # Ignore positions where both are gaps for PID
                    continue
                relevant_cols_for_pid += 1
                if c1 != '-' and c1 == c2:  # Match, and c1 is not a gap
                    matches += 1

            current_pair_pid = 0.0
            if relevant_cols_for_pid > 0:
                current_pair_pid = (matches / relevant_cols_for_pid) * 100.0

            id_matrix[i][j] = id_matrix[j][i] = current_pair_pid
            dist_matrix[i][j] = dist_matrix[j][i] = 100.0 - current_pair_pid

            total_pid_sum += current_pair_pid
            num_pairs += 1

    avg_pid = 0.0
    if num_pairs > 0:
        avg_pid = total_pid_sum / num_pairs
    elif num_seqs == 1 and aln_len > 0:  # Single sequence
        avg_pid = 100.0 if aligned_sequences[0].count('-') < aln_len else 0.0

    # 4. Sum-of-Pairs (SP) Score
    sp_score = 0.0
    if num_seqs >= 2:
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                current_pair_score = 0
                for k in range(aln_len):  # Iterate over columns
                    ci = aligned_sequences[i][k]
                    cj = aligned_sequences[j][k]
                    if ci == '-' and cj == '-':  # Gap-Gap
                        continue  # Often scored as 0 or ignored in SP calculation
                    elif ci == '-' or cj == '-':  # Gap-Residue
                        current_pair_score += gap_penalty
                    elif ci == cj:  # Match
                        current_pair_score += match_score
                    else:  # Mismatch
                        current_pair_score += mismatch_penalty
                sp_score += current_pair_score
    elif num_seqs == 1:  # SP score is typically 0 for a single sequence
        sp_score = 0.0

    return avg_pid, total_gaps_msa, fully_conserved_cols, sp_score, id_matrix, dist_matrix, num_seqs, aln_len

