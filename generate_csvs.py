"""
Generate CSV mappings between PFAM domains, proteins, and genomes.

This script processes PFAM metadata files and creates two mapping files:
1. PFAM to protein mapping
2. Protein to genome mapping
"""
import csv
import os
import sys


def generate_mappings(in_genome_master_f, in_dir, pfam2prot_f, prot2gen_f):
    """
    Generate PFAM-to-protein and protein-to-genome mapping files.

    Args:
        in_genome_master_f: Master genome file path
        in_dir: Input directory containing PFAM metadata files
        pfam2prot_f: Output PFAM-to-protein mapping file
        prot2gen_f: Output protein-to-genome mapping file
    """
    sequences = set()

    # Generate PFAM to protein mappings
    with open(pfam2prot_f, "a") as pfam2prot_file:
        filenames = next(os.walk(in_dir), (None, None, []))[2]
        for fn in filenames:
            print(f"Processing {fn}")
            with open(os.path.join(in_dir, fn), 'rt') as f:
                for line in f:
                    s = set(line.split(','))
                    sequences.update(s)
                    for prot in s:
                        pfam2prot_file.write(f"{fn},{prot}\n")

    print(f"Total unique sequences: {len(sequences)}")

    # Load genome master file
    with open(in_genome_master_f) as f:
        reader = csv.reader(f, skipinitialspace=True)
        genome_dict = dict(reader)

    # Generate protein to genome mappings
    with open(prot2gen_f, "w") as prot2gen_file:
        for protein in sequences:
            genome = genome_dict.get(protein)
            prot2gen_file.write(f"{protein},{genome}\n")

    print(f"Generated mappings: {pfam2prot_f} and {prot2gen_f}")


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python generate_csvs.py <genome_master> <input_dir> <pfam2prot_out> <prot2gen_out>")
        sys.exit(1)

    in_genome_master_f = sys.argv[1]
    in_dir = sys.argv[2]
    pfam2prot_f = sys.argv[3]
    prot2gen_f = sys.argv[4]

    generate_mappings(in_genome_master_f, in_dir, pfam2prot_f, prot2gen_f)
