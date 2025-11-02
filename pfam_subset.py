"""
Extract subset of PFAM entries from full database.

This script reads a gzipped PFAM database and extracts only the entries
matching accessions specified in an input file.
"""
import gzip


def extract_pfam_subset(pfam_file, accessions_file, output_file):
    """
    Extract PFAM entries matching specified accessions.

    Args:
        pfam_file: Path to gzipped PFAM database
        accessions_file: File containing accessions of interest (one per line)
        output_file: Output file for extracted entries
    """
    # Load accessions of interest
    with open(accessions_file, "rt") as acc_file:
        interest_acc = set(acc_file.read().splitlines())
    print(f"Loaded {len(interest_acc)} accessions of interest")

    fsm_state = "skipping"
    buffer = []
    count = 0

    with gzip.open(pfam_file, 'rt', encoding="ISO-8859-1") as f_in, \
            open(output_file, "a") as f_out:

        for line in f_in:
            if line.startswith("# STOCKHOLM"):
                fsm_state = "scanning"
                buffer = ["\n"]

            elif line.startswith("#=GF AC"):
                ac = line[7:].strip()
                if ac in interest_acc:
                    print(f"Found accession: {ac}")
                    fsm_state = "writing"
                    count += 1
                else:
                    fsm_state = "skipping"

            elif line.startswith("//") and fsm_state == "writing":
                buffer.append(line)
                print(f"Writing buffer of size {len(buffer)}")
                f_out.writelines(buffer)
                buffer = ["\n"]
                fsm_state = "reading"

            if fsm_state in ("scanning", "writing"):
                buffer.append(line)

    print(f"Extracted {count} PFAM entries")


if __name__ == '__main__':
    # Note: Update these paths as needed
    pfam_file = "/Users/arnavmeduri/Downloads/Pfam-A.full.uniprot.gz"
    accessions_file = "accessions.txt"
    output_file = "Pfam-reduced.txt"

    extract_pfam_subset(pfam_file, accessions_file, output_file)















