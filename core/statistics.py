def calculate_msa_statistics(
    aligned_sequences: list[str],
    match_score: int,
    mismatch_penalty: int,
    gap_penalty: int
) -> tuple[float, int, int, float, list[list[float]], list[list[float]], int, int]:
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
    if not aligned_sequences or not all(isinstance(s, str) for s in aligned_sequences):
        return 0.0, 0, 0, 0.0, [[]], [[]], 0, 0

    num_seqs = len(aligned_sequences)
    if num_seqs == 0:
        return 0.0, 0, 0, 0.0, [[]], [[]], 0, 0

    aln_len = len(aligned_sequences[0])
    if any(len(s) != aln_len for s in aligned_sequences):
        return 0.0, 0, 0, 0.0, \
            [[0.0] * num_seqs for _ in range(num_seqs)], \
            [[0.0] * num_seqs for _ in range(num_seqs)], \
            num_seqs, 0

    if aln_len == 0:
        empty_id_matrix = [[100.0 if i == j else 0.0 for j in range(num_seqs)] for i in range(num_seqs)]
        empty_dist_matrix = [[0.0 if i == j else 100.0 for j in range(num_seqs)] for i in range(num_seqs)]
        return 0.0, 0, 0, 0.0, empty_id_matrix, empty_dist_matrix, num_seqs, aln_len

    total_gaps_msa = _count_total_gaps(aligned_sequences)
    fully_conserved_cols = _count_fully_conserved_columns(aligned_sequences)
    avg_pid, id_matrix, dist_matrix = _calculate_pairwise_identity_and_matrices(aligned_sequences)
    sp_score = _calculate_sum_of_pairs_score(aligned_sequences, match_score, mismatch_penalty, gap_penalty)

    return (
        avg_pid,
        total_gaps_msa,
        fully_conserved_cols,
        sp_score,
        id_matrix,
        dist_matrix,
        num_seqs,
        aln_len
    )


def _count_total_gaps(aligned_sequences: list[str]) -> int:
    return sum(s.count('-') for s in aligned_sequences)


def _count_fully_conserved_columns(aligned_sequences: list[str]) -> int:
    num_seqs = len(aligned_sequences)
    aln_len = len(aligned_sequences[0])
    fully_conserved_cols = 0
    for j in range(aln_len):
        first_char_in_col = None
        is_conserved_this_col = True
        has_any_nongap_char = False
        for i in range(num_seqs):
            char_ij = aligned_sequences[i][j]
            if char_ij != '-':
                has_any_nongap_char = True
                if first_char_in_col is None:
                    first_char_in_col = char_ij
                elif first_char_in_col != char_ij:
                    is_conserved_this_col = False
                    break
        if is_conserved_this_col and has_any_nongap_char:
            fully_conserved_cols += 1
    return fully_conserved_cols


def _calculate_pairwise_identity_and_matrices(
    aligned_sequences: list[str]
) -> tuple[float, list[list[float]], list[list[float]]]:
    num_seqs = len(aligned_sequences)
    aln_len = len(aligned_sequences[0])
    total_pid_sum = 0.0
    num_pairs = 0
    id_matrix = [[0.0] * num_seqs for _ in range(num_seqs)]
    dist_matrix = [[0.0] * num_seqs for _ in range(num_seqs)]

    for i in range(num_seqs):
        if aln_len > 0 and aligned_sequences[i].count('-') < aln_len:
            id_matrix[i][i] = 100.0
        else:
            id_matrix[i][i] = 0.0
        dist_matrix[i][i] = 0.0

        for j in range(i + 1, num_seqs):
            matches = 0
            relevant_cols_for_pid = 0
            for k in range(aln_len):
                c1 = aligned_sequences[i][k]
                c2 = aligned_sequences[j][k]
                if c1 == '-' and c2 == '-':
                    continue
                relevant_cols_for_pid += 1
                if c1 != '-' and c1 == c2:
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
    elif num_seqs == 1 and aln_len > 0:
        avg_pid = 100.0 if aligned_sequences[0].count('-') < aln_len else 0.0

    return avg_pid, id_matrix, dist_matrix


def _calculate_sum_of_pairs_score(
    aligned_sequences: list[str],
    match_score: int,
    mismatch_penalty: int,
    gap_penalty: int
) -> float:
    num_seqs = len(aligned_sequences)
    aln_len = len(aligned_sequences[0])
    sp_score = 0.0
    if num_seqs >= 2:
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                current_pair_score = 0
                for k in range(aln_len):
                    ci = aligned_sequences[i][k]
                    cj = aligned_sequences[j][k]
                    if ci == '-' and cj == '-':
                        continue
                    elif ci == '-' or cj == '-':
                        current_pair_score += gap_penalty
                    elif ci == cj:
                        current_pair_score += match_score
                    else:
                        current_pair_score += mismatch_penalty
                sp_score += current_pair_score
    elif num_seqs == 1:
        sp_score = 0.0
    return sp_score