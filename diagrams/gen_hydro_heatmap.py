"""
Generate hydrophobicity heatmap for protein residue couplings.

This script analyzes amino acid couplings between proteins and generates
a heatmap showing hydrophobic vs hydrophilic interaction patterns.
"""
from collections import defaultdict
import pandas as pd


def generate_hydro_heatmap(annot_positives, protein_a='P47345', protein_b='P47346'):
    """
    Generate hydrophobicity heatmap for protein pair.

    Args:
        annot_positives: DataFrame with annotated positive couplings
        protein_a: First protein accession
        protein_b: Second protein accession

    Returns:
        Pivot table with hydrophobicity interaction counts
    """
    # Define amino acid properties
    pre = defaultdict(lambda: [])
    for aa_list, status in [
        ('EDNHQSTKR', 'hydrophilic'),
        ('ILVMFYPGW', 'hydrophobic')
    ]:
        for aa in aa_list:
            for k, v in {'aa': aa, 'type': status}.items():
                pre[k].append(v)
    aa_properties = pd.DataFrame(pre)

    # Filter and group data
    c, f = ['lft_aa', 'rgt_aa'], 'dum3'
    df = (
        annot_positives[
            (annot_positives['lft_sp_acc'] == protein_a) &
            (annot_positives['rgt_sp_acc'] == protein_b)
        ][c + [f]]
        .groupby(c, as_index=False)
        .count()
    )

    # Rename columns
    df = df.rename(columns={
        'lft_aa': protein_a,
        'rgt_aa': protein_b,
        f: 'count'
    })

    # Merge with amino acid properties
    d1 = aa_properties.rename(columns={'aa': protein_a, 'type': f'{protein_a}_type'})
    d2 = aa_properties.rename(columns={'aa': protein_b, 'type': f'{protein_b}_type'})
    df = pd.merge(left=df, how='left', on=protein_a, right=d1)
    df = pd.merge(left=df, how='left', on=protein_b, right=d2)

    # Create pivot table
    pt = df.pivot_table(
        index=f'{protein_b}_type',
        columns=f'{protein_a}_type',
        values='count',
        aggfunc=sum,
        fill_value=0
    )

    print(pt)
    return pt


# Example usage (uncomment when needed):
# if __name__ == '__main__':
#     # Load your data
#     annot_positives = pd.read_csv('path_to_your_data.csv')
#     generate_hydro_heatmap(annot_positives)
