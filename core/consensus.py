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

    for j in range(aln_len):  #columns
        column_chars = [aligned_sequences[i][j] for i in range(num_seqs)]
        counts = Counter(column_chars)

        num_gaps_in_col = counts.get('-', 0)

        #if gap frequency is above threshold, consensus is a gap
        if num_seqs > 0 and (num_gaps_in_col / num_seqs) >= gap_threshold:
            consensus_seq.append('-')
            continue

        #the most frequent non-gap character(s)
        most_common_non_gap = []
        max_freq_non_gap = 0

        #iterate through sorted characters
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
            # alphabetical tie-breaking
            most_common_non_gap.sort()
            consensus_seq.append(most_common_non_gap[0])
        else:
            # all characters in the column were gaps, but didn't meet gap_threshold
            # this shouldnt ideally not happen
            consensus_seq.append('-')

    return "".join(consensus_seq)