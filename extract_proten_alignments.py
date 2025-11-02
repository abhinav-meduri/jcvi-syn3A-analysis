"""
Extract protein alignments from PFAM Stockholm format files.

This script parses PFAM alignment files and extracts accession information.
"""
from Bio import AlignIO


def extract_protein_alignments(file_name, fmt="stockholm"):
    """
    Extract and print accession numbers from PFAM alignments.

    Args:
        file_name: Path to the alignment file
        fmt: Format of the alignment file (default: stockholm)
    """
    count = 0
    for align_it in AlignIO.parse(file_name, fmt):
        for align in align_it:
            print(align.annotations["accession"])
            count += 1
    print(f"Processed {count} alignments")


if __name__ == '__main__':
    file_name = "Pfam-reduced.txt"
    extract_protein_alignments(file_name)
