"""
Extract protein accession to organism mapping from UniProt FASTA.

This script reads a gzipped UniProt FASTA file and extracts organism (OX)
information for proteins of interest.
"""
import gzip
from Bio import SeqIO


def extract_protein_ox_mapping(uniprot_fasta, accessions_file, output_file):
    """
    Extract organism (OX) mapping for proteins of interest.

    Args:
        uniprot_fasta: Path to gzipped UniProt FASTA file
        accessions_file: File containing protein accessions of interest
        output_file: Output file for accession-OX mappings
    """
    # Load accessions of interest
    with open(accessions_file, "r") as acc_file:
        interest_acc = set(acc_file.read().splitlines())
    print(f"Loaded {len(interest_acc)} protein accessions of interest")

    count = 0
    with gzip.open(uniprot_fasta, "rt", encoding="ISO-8859-1") as handle, \
            open(output_file, "w") as f_out:

        for record in SeqIO.parse(handle, "fasta"):
            # Extract accession from ID
            acc = record.id.split('|')[1]

            if acc in interest_acc:
                # Parse description to find OX attribute
                attrs = record.description.split(" ")
                ox_attrs = [x for x in attrs if x.startswith("OX=")]

                if len(ox_attrs) == 1:
                    ox_attr_val = ox_attrs[0][3:].strip()
                    print(f"{acc}, {ox_attr_val}")
                    f_out.write(f"{acc}, {ox_attr_val}\n")
                    count += 1
                else:
                    print(f"Warning: Found {len(ox_attrs)} OX attributes for {acc}")
                    print(record.description)

    print(f"Extracted {count} protein-OX mappings")


if __name__ == '__main__':
    # Note: Update these paths as needed
    uniprot_fasta = "/Users/abhinavmeduri/Downloads/uniprot_trembl.fasta.gz"
    accessions_file = "uniq_protein_alignments.txt"
    output_file = "protein-accession-ox-mapping.txt"

    extract_protein_ox_mapping(uniprot_fasta, accessions_file, output_file)
