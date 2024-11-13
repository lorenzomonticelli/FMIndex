import sys
import os
import argparse


def parse_fasta(file_path):
    """
    The function reads a single DNA sequence from a FASTA file and adds '$' to its end.
    Args:
        file_path (str) -> Path to the FASTA file.
    Returns:
        str -> DNA sequence with a '$' terminator.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
        seq = ''.join(l.strip() for l in lines if not l.startswith('>'))
    return seq + '$'  # add terminator


def suffix_array_construction(seq):
    """
    The function builds the suffix array of a given sequence by sorting all suffixes.
    Args:
        seq (str) -> Input sequence.
    Returns:
        list -> Sorted list of starting indices of suffixes.
    """
    suffixes = sorted((seq[i:], i) for i in range(len(seq)))
    return [idx for _, idx in suffixes]


def bwt_via_sa(seq):
    """
    The function constructs the BWT of a sequence by sorting its suffixes
    and taking the last character of each sorted suffix.
    Args:
        seq (str) -> Input sequence with terminator.
    Returns:
        tuple -> BWT string and suffix array.
    """
    sa = suffix_array_construction(seq)
    return ''.join(seq[i - 1] if i != 0 else '$' for i in sa), sa


def occurrence_table(bwt):
    """
    The function builds an occurrence table for efficient rank queries in the BWT.
    Args:
        bwt (str) -> Burrows-Wheeler Transform of the sequence.
    Returns:
        dict -> Dictionary mapping each character to a list of cumulative counts across BWT.
    """
    tab = {}
    for char in set(bwt):
        tab[char] = []
        count = 0
        for sym in bwt:
            if sym == char:
                count += 1
            tab[char].append(count)
    return tab


def lf_mapping(sa, bwt):
    """
    The function creates the LF mapping array. 
    This allows traversal from BWT to the original sequence.
    Args:
        sa (list) -> Suffix array.
        bwt (str) -> Burrows-Wheeler Transform of the sequence.
    Returns:
        list -> LF mapping array.
    """
    lf = [0] * len(sa)
    first_col = sorted(bwt)
    for i in range(len(sa)):
        char = bwt[i]
        lf[i] = first_col.index(char) + bwt[:i].count(char)
    return lf


def pattern_search(pattern, bwt, occurrence_table, lf_map):
    """
    The function searches for a pattern within the FM Index using backward search.
    Args:
        pattern (str) -> The pattern to search for.
        bwt (str) -> Burrows-Wheeler Transform of the sequence.
        occurrence_table (dict) -> Occurrence table for BWT.
        lf_map (list) -> LF mapping array.
    Returns:
        list -> Indices where the pattern appears in the original sequence.
    """
    top, bottom = 0, len(bwt) - 1
    for i in range(len(pattern) - 1, -1, -1):
        char = pattern[i]
        if char in occurrence_table:
            top = occurrence_table[char][top - 1] if top > 0 else 0
            bottom = occurrence_table[char][bottom] - 1
            if top > bottom:
                return []  # pattern not found
        else:
            return []  # character not in BWT
    return list(range(top, bottom + 1))


def run_length_encode(bwt):
    """
    The function applies run-length encoding to the BWT for compression.
    Args:
        bwt (str) -> Burrows-Wheeler Transform.
    Returns:
        str -> Run-length encoded BWT.
    """
    encoded = []
    count = 1
    for i in range(1, len(bwt)):
        if bwt[i] == bwt[i - 1]:
            count += 1
        else:
            encoded.append(f"{bwt[i - 1]}{count if count > 1 else ''}")
            count = 1
    encoded.append(f"{bwt[-1]}{count if count > 1 else ''}")
    return ''.join(encoded)


def save_data(output_dir, bwt, sa, occ_table, lf, encoded_bwt=None):
    """
    The function saves the BWT, suffix array, occurrence table, LF mapping
    and encoded (if compression needed) BWT to files.
    Args:
        output_dir (str) -> Directory where files will be saved.
        bwt (str) -> Burrows-Wheeler Transform.
        sa (list) -> Suffix array.
        occ_table (dict) -> Occurrence table.
        lf (list) -> LF mapping array.
        encoded_bwt (str, optional) -> Run-length encoded BWT.
    """
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'bwt.txt'), 'w') as f:
        f.write(bwt)
    with open(os.path.join(output_dir, 'suffix_array.txt'), 'w') as f:
        f.write('\n'.join(map(str, sa)))
    with open(os.path.join(output_dir, 'occurrence_table.txt'), 'w') as f:
        for char, counts in occ_table.items():
            f.write(f"{char}: {' '.join(map(str, counts))}\n")
    with open(os.path.join(output_dir, 'lf_map.txt'), 'w') as f:
        f.write(' '.join(map(str, lf)))
    if encoded_bwt:
        with open(os.path.join(output_dir, 'bwt_rle.txt'), 'w') as f:
            f.write(encoded_bwt)


def main():
    parser = argparse.ArgumentParser(description="FM Index construction with optional pattern search and RLE.")
    parser.add_argument("fasta_file", help="Path to the FASTA file containing the DNA sequence.")
    parser.add_argument("--pattern", help="Pattern to search within the sequence.", default=None)
    parser.add_argument("--rle", action="store_true", help="Apply run-length encoding to the BWT.")

    args = parser.parse_args()

    fasta_file = args.fasta_file
    pattern = args.pattern
    apply_rle = args.rle

    # determine the directory where the FASTA file is located
    base_dir = os.path.dirname(os.path.abspath(fasta_file))
    # set the output directory
    output_dir = os.path.join(base_dir, "output")

    # process sequence
    seq = parse_fasta(fasta_file)
    
    # build FM Index components
    bwt, sa = bwt_via_sa(seq)
    occ_table = occurrence_table(bwt)
    lf_map = lf_mapping(sa, bwt)

    # apply RLE if specified
    encoded_bwt = run_length_encode(bwt) if apply_rle else None

    # save FM Index components to files
    save_data(output_dir, bwt, sa, occ_table, lf_map, encoded_bwt)
    print(f"FM Index components saved successfully in {output_dir}")

    # perform pattern search if a pattern was provided
    if pattern:
        results = pattern_search(pattern, bwt, occ_table, lf_map)
        if results:
            print(f"Pattern '{pattern}' found at positions: {results}")
        else:
            print(f"Pattern '{pattern}' not found in the sequence.")


if __name__ == "__main__":
    main()