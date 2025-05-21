from collections import Counter


def generate_consensus(aligned_sequences: list[str], gap_threshold: float = 0.51) -> str:
    """
    Generates a consensus sequence from a list of aligned sequences.

    The consensus character at each position is the most frequent non-gap character.
    If the frequency of gaps exceeds 'gap_threshold', a gap is placed in the consensus.
    Ties between non-gap characters are resolved alphabetically.

    Args:
        aligned_sequences: A list of aligned sequence strings of equal length.
        gap_threshold: The minimum fraction of gaps in a column to call a gap
                       in the consensus sequence (e.g., 0.51 means >50% gaps).

    Returns:
        The consensus sequence string.
    """
    if not aligned_sequences or not aligned_sequences[0]:
        return ""

    num_seqs = len(aligned_sequences)
    aln_len = len(aligned_sequences[0])
    consensus_seq = []

    for j in range(aln_len):  # Iterate over columns
        column_chars = [aligned_sequences[i][j] for i in range(num_seqs)]
        counts = Counter(column_chars)

        num_gaps_in_col = counts.get('-', 0)

        # Rule 1: If gap frequency is above threshold, consensus is a gap
        if num_seqs > 0 and (num_gaps_in_col / num_seqs) >= gap_threshold:
            consensus_seq.append('-')
            continue

        # Rule 2: Find the most frequent non-gap character(s)
        most_common_non_gap = []
        max_freq_non_gap = 0

        # Iterate through sorted characters to ensure deterministic tie-breaking for identical frequencies
        # (Counter.most_common() does not guarantee order for ties)
        sorted_chars = sorted(counts.keys())

        for char in sorted_chars:
            if char != '-':
                freq = counts[char]
                if freq > max_freq_non_gap:
                    max_freq_non_gap = freq
                    most_common_non_gap = [char]
                elif freq == max_freq_non_gap:
                    most_common_non_gap.append(char)

        if most_common_non_gap:
            # Alphabetical tie-breaking is implicitly handled by sorted_chars
            # if we reconstruct from that or just sort the candidates.
            # Here, `most_common_non_gap` will be built in alphabetical order if constructed carefully.
            # For simplicity, if multiple chars have max_freq_non_gap, sort them and pick the first.
            most_common_non_gap.sort()
            consensus_seq.append(most_common_non_gap[0])
        else:
            # All characters in the column were gaps, but didn't meet gap_threshold
            # This implies num_gaps_in_col < gap_threshold * num_seqs.
            # This should ideally not happen if gap_threshold is e.g. >0.5 and there's at least one non-gap,
            # or if all are gaps, it should have been caught by rule 1.
            # If it does (e.g. gap_threshold = 1.0 and all are gaps but one),
            # it might default to a gap, or N, or handle based on specific rules.
            # For now, if no non-gap character is found (e.g. all are gaps and threshold not met), append gap.
            consensus_seq.append('-')

    return "".join(consensus_seq)