"""
Generate histogram visualizations of PFAM model distributions.

This script creates visualizations showing:
1. Distribution of proteins per PFAM model per genome
2. Distribution of PFAM models per genome
"""
import pandas as pd
from matplotlib import pyplot as plt

project_dir = '.'


def load_data():
    """
    Load and merge PFAM, protein, and genome mapping data.

    Returns:
        DataFrame with merged data
    """
    # Load PFAM to protein mapping
    pfam_to_protein = pd.read_csv(
        f'{project_dir}/pfam2protein.csv',
        header=None,
        names=['pfam_id', 'protein_id']
    ).drop_duplicates()

    # Load protein to genome mapping
    protein_to_genome = pd.read_csv(
        f'{project_dir}/protein2genome.csv',
        header=None,
        names=['protein_id', 'genome_id']
    ).drop_duplicates()

    # Merge tables
    big_table = pd.merge(
        left=pfam_to_protein,
        right=protein_to_genome,
        on='protein_id'
    )

    big_table.to_csv('big_table.csv')
    print("Created big_table.csv")
    return big_table


big_table = load_data()

# Count proteins per PFAM model per genome
protein_per_genome_count = (
    big_table.groupby(['genome_id', 'pfam_id'], as_index=False)
    .count()
    .rename(columns={'protein_id': 'protein_count'})
)

# Plot distribution of proteins per model+genome
plt.figure(figsize=(20, 10))
max_count = max(protein_per_genome_count['protein_count'])
print(f"Max protein_per_genome_count: {max_count}")
plt.hist(protein_per_genome_count['protein_count'], bins=range(0, 200))
plt.xlabel('Proteins per model + genome')
plt.ylabel('Frequency')
plt.title('Distribution of Proteins per PFAM Model per Genome')
plt.grid(True, alpha=0.3)
plt.show()

# Plot distribution of number of models per genome
c, f = ['genome_id'], 'pfam_id'
model_per_genome_count = (
    protein_per_genome_count[c + [f]]
    .groupby(c, as_index=False)
    .count()
    .rename(columns={f: 'model_count'})
)

plt.figure(figsize=(20, 10))
mx = model_per_genome_count['model_count']
print(f"Number of genomes: {len(mx)}")
print(f"Max models per genome: {int(max(mx))}")

# Filter dataset for genomes with >= cutoff models
cutoff = 10
df = model_per_genome_count[model_per_genome_count['model_count'] >= cutoff]
plt.hist(df['model_count'], bins=100)
plt.xlabel('PFAM models per genome')
plt.ylabel('Frequency')
plt.title(f'Distribution of PFAM Models per Genome (>= {cutoff} models)')
plt.grid(True, alpha=0.3)
plt.show()
