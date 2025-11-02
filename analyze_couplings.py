"""
Analyze protein coupling scores from PFAM alignments.

This script processes coupling score files and filters them based on a cutoff threshold.
"""
import sys
import os
import pandas as pd


def analyze_couplings(filepath, output_file_name):
    """
    Process coupling scores and filter based on threshold.

    Args:
        filepath: Path to the coupling scores file
        output_file_name: Path to the output CSV file
    """
    print(f"processing file: {filepath}")

    d, f = os.path.split(filepath)
    pfams = f.split('.')[0].strip()
    lft_pfam_id = pfams.split('-')[0].strip()
    rgt_pfam_id = pfams.split('-')[1].strip()

    couplings = pd.read_csv(
        filepath,
        sep=" ",
        header=None,
        names=["pos_lo", "dum1", "dum2", "pos_hi", "dum3", "coupling_score"]
    )[['pos_lo', 'pos_hi', 'coupling_score']]

    # Add columns for the left and right PFAM IDs
    couplings['pfam_lft'] = lft_pfam_id
    couplings['pfam_rgt'] = rgt_pfam_id

    # Filter couplings based on cutoff threshold
    # Higher coupling scores indicate stronger evidence of correlation
    cutoff = 0.0001
    filtered = couplings[couplings['coupling_score'] >= cutoff]
    filtered.to_csv(output_file_name, header=False, index=False, mode='a')
    print(f"wrote into: {output_file_name}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python analyze_couplings.py <input_file> <output_file>")
        sys.exit(1)

    filepath = sys.argv[1]
    output_file_name = sys.argv[2]
    analyze_couplings(filepath, output_file_name)
