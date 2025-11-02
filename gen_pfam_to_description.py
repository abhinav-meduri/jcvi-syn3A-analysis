"""
Generate PFAM ID to description mapping from UniProt database.

This script parses a gzipped UniProt database in Stockholm format and extracts
PFAM accessions along with their descriptions.
"""
import gzip

import click


@click.command()
@click.option('--uniprot-db', required=True, help='Path to gzipped UniProt database file')
@click.option('--out-file', required=True, help='Output CSV file path')
def scan_db(uniprot_db, out_file):
    """
    Scan UniProt database and extract PFAM ID to description mappings.

    Args:
        uniprot_db: Path to gzipped UniProt database file
        out_file: Output CSV file path
    """
    with open(out_file, "w") as f_out:
        with gzip.open(uniprot_db, 'rt', encoding='ISO-8859-1') as f:
            pfam_id = desc = ""

            for line in f:
                if line.startswith("# STOCKHOLM"):
                    # Reset accession and description strings
                    pfam_id = desc = ""

                elif line.startswith("#=GF AC"):
                    if pfam_id:
                        print(f"Warning: pfam_id {pfam_id} already set")
                    pfam_id = line[7:].strip().split('.')[0]

                elif line.startswith("#=GF DE"):
                    if desc:
                        print(f"Warning: description {desc} already set")
                    desc = line[7:].strip()

                elif line.startswith("//"):
                    if pfam_id and desc:
                        f_out.write(f"{pfam_id},{desc}\n")
                    else:
                        print(f"Warning: Missing pfam_id or description")


if __name__ == '__main__':
    scan_db()
