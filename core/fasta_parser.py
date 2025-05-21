def parse_fasta(file_content: str, max_header_length: int = 10) -> list[tuple[str, str]]:
    """
    Parses a string containing FASTA formatted sequences.

    Args:
        file_content: The string content of the FASTA file.
        max_header_length: The maximum allowed length for sequence headers.
                           Headers longer than this will be truncated.

    Returns:
        A list of tuples, where each tuple is (header, sequence_string).
        Headers are truncated if necessary and made unique if duplicates are found.
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
                # Add previous sequence (header already processed for length)
                sequences.append((current_header, "".join(current_sequence_parts).upper()))

            # Extract header, strip whitespace, and then truncate
            raw_header = line[1:].strip()
            current_header = raw_header[:max_header_length]  # Apply length limit
            current_sequence_parts = []
        elif current_header is not None:  # Sequence line
            current_sequence_parts.append(line.replace(" ", "").replace("\t", ""))  # Remove spaces/tabs within sequence

    # Add the last sequence
    if current_header is not None and current_sequence_parts:
        # Header (current_header) was already processed for length when it was read
        sequences.append((current_header, "".join(current_sequence_parts).upper()))

    # Ensure unique headers (operating on potentially truncated headers)
    final_sequences = []
    used_headers = set()
    for idx, (header, seq) in enumerate(sequences):
        # If header became empty after stripping/truncation, provide a default
        original_header_processed = header if header else f"UnnamedSeq{idx + 1}"
        # Ensure the default name also respects max_header_length, though unlikely to be an issue for "UnnamedSeq"
        original_header_processed = original_header_processed[:max_header_length]

        current_name_candidate = original_header_processed
        count = 1
        # Check for uniqueness among (potentially truncated) headers
        while current_name_candidate in used_headers:
            # If truncated name conflicts, append suffix to the *truncated* base
            base_for_suffix = original_header_processed
            suffix = f"_{count}"

            # Ensure the new name with suffix doesn't exceed max_header_length
            # If it does, truncate the base part further to make space for the suffix
            if len(base_for_suffix) + len(suffix) > max_header_length:
                base_for_suffix = base_for_suffix[:max_header_length - len(suffix)]

            current_name_candidate = f"{base_for_suffix}{suffix}"
            count += 1

        used_headers.add(current_name_candidate)
        final_sequences.append((current_name_candidate, seq))

    return final_sequences

