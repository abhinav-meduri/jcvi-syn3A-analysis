"""
Generate PyMOL scripts for visualizing protein-protein interaction surfaces.

This script creates PyMOL script files (.pml) that highlight residues involved
in protein-protein interactions with consistent color coding.
"""
import click

# Color palette for residue highlighting
COLORS = [
    "red", "green", "blue", "yellow", "magenta", "cyan", "orange", "salmon",
    "palegreen", "lightblue", "wheat", "lightpink", "lightteal", "lightorange",
    "brown", "limon", "purpleblue", "paleyellow", "hotpink", "palecyan", "olive"
]

# PyMOL script templates
TEMPLATE_PROTEIN_1 = """bg_color white
load L_SP.pdb
cmd.color_deep("gray90", 'L_SP', 0)
"""

TEMPLATE_PROTEIN_2 = """
load R_SP.pdb
cmd.color_deep("gray40", 'R_SP', 0)
"""

TEMPLATE_RESIDUE = """
select KEEPER, SP and resi POS
color COLOR, KEEPER
show spheres, KEEPER
label name CA+C1*+C1' and byres(KEEPER), "%s-%s"% (resn,resi)
"""


@click.command()
@click.option('--hits_file', default='./surfaces.csv',
              help='CSV file with protein interaction surfaces')
def scan_args(hits_file):
    """
    Generate PyMOL scripts from surface interaction file.

    Args:
        hits_file: CSV file with columns: l_pos,l_sp,r_pos,r_sp
    """
    with open(hits_file, 'r') as f:
        for line in f.readlines():
            l_pos, l_sp, r_pos, r_sp = line.strip().split(',')
            l_pos_list = l_pos.split('+')
            r_pos_list = r_pos.split('+')

            assert len(l_pos_list) == len(r_pos_list), \
                f"Mismatch in position lists: {len(l_pos_list)} vs {len(r_pos_list)}"

            # Track colors for consistent pairing
            l_dict = {}
            color_idx = 0

            output_file = f'{l_sp}-{r_sp}-c.pml'
            with open(output_file, 'w') as f_out:
                # Setup left protein
                f_out.write(TEMPLATE_PROTEIN_1.replace("L_SP", l_sp))

                # Highlight left protein residues
                for acid in l_pos_list:
                    if acid in l_dict:
                        color = l_dict[acid]
                    else:
                        color = COLORS[color_idx]
                        color_idx += 1
                        l_dict[acid] = color

                    residue_cmd = (
                        TEMPLATE_RESIDUE
                        .replace("SP", l_sp)
                        .replace("KEEPER", f"l_{acid}")
                        .replace("COLOR", color)
                        .replace("POS", acid)
                    )
                    f_out.write(residue_cmd)

                # Setup right protein
                f_out.write(TEMPLATE_PROTEIN_2.replace("R_SP", r_sp))

                # Highlight right protein residues with matching colors
                for x, acid in enumerate(r_pos_list):
                    color = l_dict[l_pos_list[x]]
                    residue_cmd = (
                        TEMPLATE_RESIDUE
                        .replace("SP", r_sp)
                        .replace("KEEPER", f"r_{acid}")
                        .replace("COLOR", color)
                        .replace("POS", acid)
                    )
                    f_out.write(residue_cmd)

            print(f"Generated {output_file}")


if __name__ == '__main__':
    scan_args()
