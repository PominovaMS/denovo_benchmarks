"""Script to convert predicted labels from the original algorithm 
output format to the common data format."""

import argparse
import re
from pyteomics.mztab import MzTab


REPLACEMENTS = [
    ("C+57.021", "C")  # C is written without Carbamidomethyl modification
]

PTM_PATTERN = r"^([0-9.+-]+)([A-Z])"


def transform_match(match: re.Match) -> str:
    """
    Transform representation of amino acids substring matching
    the PTM pattern.

    Parameters
    ----------
    match : re.Match
        Substring matching the PTM pattern.

    Returns
    -------
    transformed_match : str
        Transformed PTM pattern representation.
    """
    ptm, aa = match.group(1), match.group(2)
    return aa + ptm


def convert_to_output_format(sequence: str) -> str:
    """
    Convert peptide sequence to the common output data format.

    Parameters
    ----------
    sequence : str
        Peptide sequence in the original algorithm output format.

    Returns
    -------
    transformed_sequence : str
        Peptide sequence in the common output data format.
    """

    # direct (token-to-token) replacements
    for repl_args in REPLACEMENTS:
        sequence = sequence.replace(*repl_args)

    # transformation of PTM notation
    # move N-terminus modifications BEYOND 1st AA
    sequence = re.sub(PTM_PATTERN, transform_match, sequence)

    return sequence


parser = argparse.ArgumentParser()
parser.add_argument(
    "output_path", help="The path to the algorithm predictions file."
)
args = parser.parse_args()

# read predictions from output file
output_data = MzTab(args.output_path)

output_data = output_data.spectrum_match_table
output_data["sequence"] = output_data["sequence"].apply(
    convert_to_output_format
)

# save processed predictions to the same file
output_data.to_csv("outputs.csv", index=False)
