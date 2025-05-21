def parse_fasta(file_content: str, max_header_length: int = 10) -> list[tuple[str, str]]:
    """
    Parses a string containing FASTA formatted sequences.

    Args:
        file_content: The string content of the FASTA file.
        max_header_length: The maximum allowed length for sequence headers.
                           Headers longer than this will be truncated.

    Returns:
        A list of tuples, where each tuple is (a header, sequence_string).
        Headers are truncated if necessary and made unique if duplicates are found.
    """
    sequences = []
    current_header = None
    current_sequence_parts = []

    for line in file_content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None and current_sequence_parts:
                sequences.append((current_header, "".join(current_sequence_parts).upper()))

            raw_header = line[1:].strip()
            current_header = raw_header[:max_header_length]
            current_sequence_parts = []
        elif current_header is not None:
            current_sequence_parts.append(line.replace(" ", "").replace("\t", ""))

    if current_header is not None and current_sequence_parts:
        sequences.append((current_header, "".join(current_sequence_parts).upper()))

    final_sequences = []
    used_headers = set()
    for idx, (header, seq) in enumerate(sequences):
        original_header_processed = header if header else f"UnnamedSeq{idx + 1}"
        original_header_processed = original_header_processed[:max_header_length]

        current_name_candidate = original_header_processed
        count = 1
        while current_name_candidate in used_headers:
            base_for_suffix = original_header_processed
            suffix = f"_{count}"

            if len(base_for_suffix) + len(suffix) > max_header_length:
                base_for_suffix = base_for_suffix[:max_header_length - len(suffix)]

            current_name_candidate = f"{base_for_suffix}{suffix}"
            count += 1

        used_headers.add(current_name_candidate)
        final_sequences.append((current_name_candidate, seq))

    return final_sequences

