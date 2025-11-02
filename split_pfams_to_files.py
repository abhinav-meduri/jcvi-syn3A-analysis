"""
Split PFAM database into individual files per accession.

This script reads a gzipped PFAM database and splits it into separate files,
one for each PFAM accession specified in the input file.

Usage:
    python split_pfams_to_files.py Pfam-A.full.uniprot.gz JCVI-Syn3A-pfam-accessions <target-dir>
"""
import sys
import os
import gzip


def split_pfam(file_name, acc_f, dir_path):
    """
    Split PFAM database into individual files.

    Args:
        file_name: Path to gzipped PFAM database
        acc_f: File containing PFAM accessions of interest
        dir_path: Output directory for individual PFAM files
    """
    # Load accessions
    with open(acc_f, "rt") as acc_file:
        interest_acc = acc_file.read().splitlines()

    count = len(interest_acc)
    print(f"Total number of PFAM accessions to process: {count}")

    buffer = []
    fsm_state = "undefined"

    with gzip.open(file_name, "rt", encoding="ISO-8859-1") as f:
        for line in f:
            if line.startswith("# STOCKHOLM"):
                fsm_state = "scanning"
                buffer = [line]
                sequences = set()

            elif line.startswith("#=GF AC"):
                ac = line[7:].strip().split(".", 1)[0]
                if ac in interest_acc:
                    assert fsm_state == "scanning", f"Unexpected state: {fsm_state}"
                    fsm_state = "writing"
                    buffer.append(line)
                    print(f"Scanning PFAM: {ac}")
                else:
                    fsm_state = "skipping"

            elif line[0] != "#" and not line.startswith("//") and fsm_state == "writing":
                buffer.append(line)
                # Sequence line format: seq-name sequence
                parts = [x.strip() for x in line.split(" ", 1)]
                if len(parts) != 2:
                    raise ValueError(
                        f"Could not split line into identifier and sequence\n{line}"
                    )
                seq_id = parts[0].split('/')[0].split('.')[0]
                sequences.add(seq_id)

            elif line.startswith("//") and fsm_state == "writing":
                buffer.append(line)

                # Write PFAM file
                out_file_name = os.path.join(dir_path, ac)
                with open(out_file_name, "w") as out_file:
                    out_file.writelines(buffer)

                # Write metadata file
                meta_dir = os.path.join(dir_path, "meta")
                os.makedirs(meta_dir, exist_ok=True)
                with open(os.path.join(meta_dir, ac), "w") as out_file_meta:
                    out_file_meta.write(','.join(sequences))

                # Reset state
                buffer = []
                sequences = set()
                fsm_state = "undefined"
                count -= 1
                print(f"Wrote PFAM {ac}, {count} more to go")

                if count <= 0:
                    print("Done!")
                    return

            elif fsm_state == "found_ac":
                buffer.append(line)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python split_pfams_to_files.py <pfam_db.gz> <accessions_file> <output_dir>")
        sys.exit(1)

    input_pfam_file = sys.argv[1]
    input_acc_file = sys.argv[2]
    out_dir = sys.argv[3]

    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)

    split_pfam(input_pfam_file, input_acc_file, out_dir)
