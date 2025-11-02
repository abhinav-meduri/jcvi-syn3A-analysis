"""
Generate FASTA file from GenBank format genome file.

This script parses a GenBank file and extracts all CDS (coding sequence) features,
generating a FASTA file with protein sequences and metadata.
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_cds_features(file_name):
    """
    Extract CDS features from GenBank file and generate FASTA.

    Args:
        file_name: Path to GenBank file
    """
    # Stores all the CDS entries
    all_entries = []

    gb_record = SeqIO.read(open(file_name, "r"), "genbank")
    gb_feature = gb_record.features
    accession = (gb_record.annotations['accessions'][0] + '.' +
                 str(gb_record.annotations['sequence_version']))

    chromosome = "not-defined"
    strain = "not-defined"
    organism = "not-defined"

    # Extract organism, strain and chromosome values
    for seq_feature in gb_feature:
        if seq_feature.type == "source":
            organism = str(seq_feature.qualifiers['organism'][0]).rsplit(" ", 1)[0]
            if 'strain' in seq_feature.qualifiers:
                strain = seq_feature.qualifiers['strain'][0]
            if 'chromosome' in seq_feature.qualifiers:
                chromosome = seq_feature.qualifiers['chromosome'][0]

    # Process all CDS features
    for seq_feature in gb_feature:
        if seq_feature.type == "CDS":
            if 'translation' in seq_feature.qualifiers:
                if seq_feature.qualifiers['translation'][0]:
                    # Extract various qualifiers with defaults
                    gene = seq_feature.qualifiers.get('gene', ['not-defined'])[0]
                    product = seq_feature.qualifiers.get('product', ['not-defined'])[0]
                    protein_id = seq_feature.qualifiers.get('protein_id', ['not-defined'])[0]
                    old_locus_tag = seq_feature.qualifiers.get('old_locus_tag', ['not-defined'])[0]
                    locus_tag = seq_feature.qualifiers.get('locus_tag', ['not-defined'])[0]

                    # Determine strand complement
                    if seq_feature.location.strand == 1:
                        complement = 'no-complement'
                    elif seq_feature.location.strand == -1:
                        complement = str(seq_feature.location).strip('[]').split('(-)')[0]
                    else:
                        complement = 'not-defined'

                    # Create SeqRecord with metadata
                    seq_record = SeqRecord(
                        Seq(seq_feature.qualifiers['translation'][0]),
                        id=locus_tag,
                        description=(
                            f"|{old_locus_tag}|{gene}|{product}|{complement}|{protein_id}|"
                            f"Chromosome-{chromosome}|{accession}|{organism}|{strain}"
                        )
                    )
                    all_entries.append(seq_record)

    # Write FASTA file
    output_file = f'{file_name[:-3]}.fasta'
    SeqIO.write(all_entries, output_file, 'fasta')
    print(f"Generated {output_file} with {len(all_entries)} CDS entries")


if __name__ == '__main__':
    file_name = 'JCVI-Syn3A.gb'
    extract_cds_features(file_name)
