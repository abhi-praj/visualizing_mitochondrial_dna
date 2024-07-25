def validate_file(fasta_file: str) -> bool:
    valid_dna_characters = set("ATCGRYMKSWHBVDN")

    ambiguity_map = {
        'R': {'A', 'G'},  # R = A or G
        'Y': {'C', 'T'},  # Y = C or T
        'M': {'A', 'C'},  # M = A or C
        'K': {'G', 'T'},  # K = G or T
        'S': {'C', 'G'},  # S = C or G
        'W': {'A', 'T'},  # W = A or T
        'H': {'A', 'C', 'T'},  # H = A or C or T
        'B': {'C', 'G', 'T'},  # B = C or G or T
        'V': {'A', 'C', 'G'},  # V = A or C or G
        'D': {'A', 'G', 'T'},  # D = A or G or T
        'N': {'A', 'T', 'C', 'G'}  # N = any base
    }

    is_header = True

    with open(fasta_file, 'r') as f:
        sequence_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                is_header = True
                if sequence_lines:
                    sequence = ''.join(sequence_lines)
                    if any(base not in valid_dna_characters for base in
                           sequence.upper()):
                        return False
                    sequence_lines = []
            else:
                if is_header:
                    is_header = False
                for code, bases in ambiguity_map.items():
                    for base in bases:
                        line = line.replace(code, base)
                sequence_lines.append(line)

        if sequence_lines:
            sequence = ''.join(sequence_lines)
            if any(base not in valid_dna_characters for base in
                   sequence.upper()):
                return False

    return True
