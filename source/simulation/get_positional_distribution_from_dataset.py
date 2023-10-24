'''
The purpose of this script is to retrieve the positional amino acid distribution from a given dataset (e.g., Mason),
and write it to an output file. This output file can then be used by the script simulate_synthetic_data.py in order
to simulate more data according to the same positional frequency distribution.
'''

import sys
import argparse
from source.utility import *

def export_positional_distribution(np_sequences, output_file):
    data = {"position": [],
            "amino_acid": [],
            "relative_frequency": []}

    for pos in range(len(np_sequences[0])):

        for amino_acid in ALPHABET:
            data["position"].append(pos)
            data["amino_acid"].append(amino_acid)
            data["relative_frequency"].append(sum(np_sequences[:, pos] == amino_acid))

    df = pd.DataFrame(data)
    df.to_csv(output_file, sep="\t", index=False)


def parse_commandline_arguments(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("--data_file", default="../../data/mason/train_val_test/mason_all.csv", type=str)
    parser.add_argument("--output_file", default="../../data/simulated_data/positional_distribution_mason.tsv", type=str)
    parser.add_argument("--log_file", default="../../data/simulated_data/log_positional_distribution_mason.tsv", type=str)

    parsed_args = parser.parse_args(args)
    return parsed_args


if __name__ == "__main__":
    args = parse_commandline_arguments(sys.argv[1:])

    log_arguments(args.log_file, args)

    np_sequences, _ = read_data_file(args.data_file, use_np_repr=True)

    export_positional_distribution(np_sequences, args.output_file)

