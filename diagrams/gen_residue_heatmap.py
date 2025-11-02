"""
Generate residue-level heatmap for protein couplings.

This script visualizes protein-protein interaction couplings at the residue level
using heatmaps based on z-score significance.
"""
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def generate_residue_heatmap(coupling_file, project_dir='.', zscore_cutoff=4):
    """
    Generate residue-level heatmap for protein couplings.

    Args:
        coupling_file: Path to coupling annotations CSV file
        project_dir: Project directory path
        zscore_cutoff: Z-score threshold for filtering

    Returns:
        Pivot table and displays heatmap
    """
    # Load coupling data
    hit1_subset = pd.read_csv(f'{project_dir}/{coupling_file}')

    # Filter significant inter-protein couplings
    hit1_subset = hit1_subset[
        (hit1_subset['zscore'] >= zscore_cutoff) &
        (hit1_subset['lft_sp_id'] != hit1_subset['rgt_sp_id'])
    ][[
        'zscore', 'lft_sp_id', 'rgt_sp_id',
        'lft_aa', 'lft_sp_pos', 'rgt_aa', 'rgt_sp_pos'
    ]]

    print(f"Unique left positions: {len(hit1_subset['lft_sp_pos'].unique())}")
    print(f"Unique right positions: {len(hit1_subset['rgt_sp_pos'].unique())}")
    print(f"Total interactions: {hit1_subset.shape[0]}")

    # Create pivot table
    pt = hit1_subset.pivot_table(
        index=['lft_sp_pos'],
        columns=['rgt_sp_pos'],
        values='zscore',
        aggfunc=len
    )

    # Generate heatmap
    plt.figure(figsize=(20, 10))
    sns.heatmap(pt, cmap='viridis')
    plt.xlabel('Right protein position')
    plt.ylabel('Left protein position')
    plt.title(f'Residue-level Coupling Heatmap (z-score >= {zscore_cutoff})')
    plt.tight_layout()
    plt.show()

    return pt


# Example usage:
# if __name__ == '__main__':
#     generate_residue_heatmap('PF01425-PF02934.couplings.norm.annot.csv')
