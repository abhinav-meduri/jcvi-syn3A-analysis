"""
Filter PFAM results based on log10 cutoff thresholds.

This script processes PFAM search results and filters them based on e-value cutoffs,
generating plots to visualize the relationship between cutoff values and number of hits.
"""
import pandas as pd
from matplotlib import pyplot as plt

project_dir = '.'

cnames = ''.join([
    '', 'target_name,prot_accession,tlen,pfam_name,pfam_id_versioned,qlen,full_e_value,',
    'full_score,full_bias,domain_n,domain_of,domain_c_evalue,',
    'domain_i_evalue,domain_score,domain_bias,',
    'hmm_from,coord_to_1,ali_from,coord_to_2,env_from,coord_to_3,acc,', 'target_description',
]).split(',')

# Load PFAM results
gsheet = pd.read_csv(
    f'{project_dir}/PFAM-results.xlsx - original.csv',
    skiprows=2,  # Skip first two rows - ignore current headers
    comment='#',  # Ignore lines starting with # - comments at bottom of file
    header=None,
    names=cnames  # New column names
)

# Shape of a dataframe is in form rows, columns
print(f"Data shape: {gsheet.shape}")

# Display key columns
print(gsheet[[
    'target_name', 'tlen', 'pfam_name', 'pfam_id_versioned',
    'qlen', 'full_e_value', 'target_description'
]])

# Strip version number from pfam_id_versioned and create new column
# apply runs a function on each row (axis=1 for rows, axis=0 for columns)
gsheet['pfam_id'] = gsheet.apply(
    lambda row: row['pfam_id_versioned'].split('.')[0],
    axis=1
)

# Filter on e-value
gsheet_filtered = gsheet[gsheet['full_e_value'] < 1e-02]
print(f"Filtered shape: {gsheet_filtered.shape}")
gsheet_filtered[['pfam_id']].to_csv("/tmp/a.csv", index=False)

print(f"Number of unique PFAM IDs: {len(gsheet_filtered['pfam_id'].unique())}")

c, f = ['pfam_id'], 'target_name'
cutoff = 1e-02

# Count proteins per PFAM model
myco_protein_counts = (
    gsheet[gsheet['full_e_value'] <= cutoff][c + [f]]
    .groupby(c, as_index=False)
    .count()
    .rename(columns={f: 'protein_count'})
)

# How many models match only once in Mycoplasma proteome?
single_matches = myco_protein_counts[myco_protein_counts['protein_count'] == 1].shape[0]
print(f"Number of models matched only once in Mycoplasma proteome: {single_matches}")

# Visualize effect of different cutoffs
plt.figure(figsize=(20, 4))
x, y = [], []
for exponent in range(-2, -90, -1):
    x.append(exponent)
    df = gsheet[gsheet['full_e_value'] <= 10 ** exponent]
    y.append(len(df['pfam_id'].unique()))
plt.plot(x, y)
plt.xlabel('log10 cutoff')
plt.ylabel('Number of PFAM models with >=1 hit')
plt.title('PFAM Models vs E-value Cutoff')
plt.grid(True, alpha=0.3)
plt.show()
