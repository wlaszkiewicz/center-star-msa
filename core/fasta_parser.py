def parse_fasta(file_content: str) -> list[tuple[str, str]]:
    """
    Parses a string containing FASTA formatted sequences.

    Args:
        file_content: The string content of the FASTA file.

    Returns:
        A list of tuples, where each tuple is (header, sequence_string).
        Headers are made unique if duplicates are found.
    """
    sequences = []
    current_header = None
    current_sequence_parts = []

    for line in file_content.splitlines():
        line = line.strip()
        if not line:  # Skip empty lines
            continue
        if line.startswith(">"):
            if current_header is not None and current_sequence_parts:
                sequences.append((current_header, "".join(current_sequence_parts).upper()))
            current_header = line[1:].strip()  # Remove ">" and strip whitespace
            current_sequence_parts = []
        elif current_header is not None:  # Sequence line
            current_sequence_parts.append(line.replace(" ", "").replace("\t", ""))  # Remove spaces/tabs within sequence

    # Add the last sequence
    if current_header is not None and current_sequence_parts:
        sequences.append((current_header, "".join(current_sequence_parts).upper()))

    # Ensure unique headers
    final_sequences = []
    used_headers = set()
    for idx, (header, seq) in enumerate(sequences):
        original_header = header if header else f"UnnamedSeq{idx + 1}"
        current_name_candidate = original_header
        count = 1
        while current_name_candidate in used_headers:
            current_name_candidate = f"{original_header}_{count}"
            count += 1
        used_headers.add(current_name_candidate)
        final_sequences.append((current_name_candidate, seq))

    return final_sequences