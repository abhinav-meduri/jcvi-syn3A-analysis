"""
Glue together PFAM alignments to create concatenated alignments.

This script combines alignments from multiple PFAM domains that co-occur in
the same genomes, creating glued alignment files for evolutionary coupling analysis.
"""
import os
import sys
import time
import functools

import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def timer(func):
    """
    Decorator to time function execution.

    Args:
        func: Function to time

    Returns:
        Wrapped function with timing
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        print(f"Finished {func.__name__!r} in {round(run_time, 3)} secs")
        return value
    return wrapper


def parse_protein_name(alignid):
    """
    Extract protein name from alignment ID.

    Args:
        alignid: Alignment identifier

    Returns:
        Protein name without version
    """
    return alignid.split('.')[0]


@timer
def stockholm_to_dataframe(pfam_id, protein_to_species):
    """
    Convert Stockholm format PFAM alignment to filtered dataframe.

    Args:
        pfam_id: PFAM accession identifier
        protein_to_species: DataFrame mapping proteins to genome IDs

    Returns:
        Filtered DataFrame with alignments unique to each genome
    """
    mult = AlignIO.read(f'{stockholm_files_dir}/{pfam_id}', 'stockholm')
    pre = {k: [] for k in ('protein_id', 'protein_name', 'alignment')}

    for indiv in mult:
        # Each alignment in PFAM file is derived from one protein
        # Parse protein ID from alignment identifier
        protein_id = parse_protein_name(indiv.id)

        # Store protein information
        for k, v in {
            'protein_id': protein_id,
            'protein_name': indiv.id,
            'alignment': str(indiv.seq)
        }.items():
            pre[k].append(v)

    # Create initial dataframe
    df0 = pd.DataFrame(pre)
    df0['pfam_id'] = pfam_id

    # Merge with genome information (inner join)
    full_df = pd.merge(
        left=df0,
        right=protein_to_species[['protein_id', 'genome_id']],
        on=['protein_id']
    )

    # Compute counts by genome
    c, f = ['genome_id'], 'protein_id'
    genome_counts = (
        full_df[c + [f]]
        .groupby(c, as_index=False)
        .count()
        .rename(columns={f: 'count_for_genome'})
    )

    # Filter for proteins unique to each genome
    filtered_df = pd.merge(
        left=genome_counts[genome_counts['count_for_genome'] == 1][['genome_id']],
        right=full_df,
        on='genome_id'
    )

    return filtered_df


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python glue_alignments.py <stockholm_dir> <out_dir> <protein_to_species> <pfam_id>")
        sys.exit(1)

    project_dir = '.'
    stockholm_files_dir = sys.argv[1]
    out_dir = sys.argv[2]
    protein_to_species_f = sys.argv[3]
    pfam_id_a = sys.argv[4]

    # Load protein to species mapping
    protein_to_species = pd.read_csv(
        protein_to_species_f,
        header=None,
        names=['protein_id', 'genome_id']
    )

    # Load PFAM accessions
    pfam_set = set(line.strip() for line in open(f'{stockholm_files_dir}/accessions.csv'))
    if pfam_id_a in pfam_set:
        pfam_set.remove(pfam_id_a)
        print(f"Processing {pfam_id_a} against {len(pfam_set)} other PFAMs")
    else:
        print(f"PFAM {pfam_id_a} not in set!")
        sys.exit(1)

    # Get alignments for the primary PFAM
    align_df_a = stockholm_to_dataframe(pfam_id_a, protein_to_species)

    # Process each pair
    for pfam_id_b in pfam_set:
        align_df_b = stockholm_to_dataframe(pfam_id_b, protein_to_species)

        # Merge alignments by genome
        p_df = pd.merge(
            left=align_df_a,
            right=align_df_b,
            on=['genome_id'],
            suffixes=['.lft', '.rgt']
        )

        # Check for Mycoplasma genitalium (genome_id 243273)
        mycge = p_df[p_df['genome_id'] == 243273]
        if mycge.shape[0] == 0:
            print("Did not find proteins related to 243273, skipping")
            continue

        # Save Mycoplasma-specific results
        os.makedirs(f'{out_dir}/glued-{pfam_id_a}', exist_ok=True)
        mycge.to_csv(f'{out_dir}/glued-{pfam_id_a}/{pfam_id_a}-{pfam_id_b}.mycge.csv')

        # Create glued sequence records
        seq_records = []
        for ridx, row in p_df.iterrows():
            seq_records.append(
                SeqRecord(
                    Seq(row['alignment.lft'] + row['alignment.rgt']),
                    id=row['protein_name.lft'] + '+' + row['protein_name.rgt']
                )
            )

        # Write glued alignments if sufficient sequences
        num_seqs = len(seq_records)
        if num_seqs > 8192:
            print(f"\tWriting {num_seqs} sequences")
            with open(f'{out_dir}/glued-{pfam_id_a}/{pfam_id_a}-{pfam_id_b}.faln', 'w') as wrt:
                SeqIO.write(seq_records, wrt, 'fasta')
        else:
            print(f"Glued alignments [{pfam_id_a} {pfam_id_b}] have only {num_seqs} sequences, skipping")
