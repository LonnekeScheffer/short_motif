import sys
import random
import argparse
from source.utility import *


def get_amino_acids_from_positional_distribution(indices, positional_distribution):
    if positional_distribution is None:
        return random.choices(ALPHABET, k=len(indices))
    else:
        return [random.choices(population=positional_distribution[index]["amino_acids"],
                               weights=positional_distribution[index]["relative_frequencies"],
                               k=1)[0] for index in indices]

def generate_random_motif(motif_size, sequence_length=SEQUENCE_LENGTH, positional_distribution=None):
    motif_indices = sorted(random.sample(range(sequence_length), k=motif_size))

    amino_acids = get_amino_acids_from_positional_distribution(motif_indices, positional_distribution)

    return (tuple(motif_indices), tuple(amino_acids))

def generate_random_motifs(n_motifs, motif_size, sequence_length=SEQUENCE_LENGTH, positional_distribution=None):
    return [generate_random_motif(motif_size, sequence_length, positional_distribution) for i in range(n_motifs)]

def generate_random_positive_sequence(sequence_length, implant_motif, positional_distribution=None):
    non_motif_indices = [i for i in range(sequence_length) if i not in implant_motif[0]]
    sequence = get_amino_acids_from_positional_distribution(non_motif_indices, positional_distribution)

    for index, amino_acid in zip(implant_motif[0], implant_motif[1]):
        sequence.insert(index, amino_acid)

    return "".join(sequence)

def test_contains_motif(sequence, motif):
    for index, amino_acid in zip(motif[0], motif[1]):
        if sequence[index] != amino_acid:
            return False

    return True


def test_contains_any_motif(sequence, motifs):
    for motif in motifs:
        if test_contains_motif(sequence, motif):
            return True

    return False

def generate_random_negative_sequence(sequence_length, motifs, positional_distribution=None):
    sequence = "".join(get_amino_acids_from_positional_distribution(range(sequence_length), positional_distribution))

    while test_contains_any_motif(sequence, motifs):
        sequence = "".join(get_amino_acids_from_positional_distribution(range(sequence_length), positional_distribution))

    return sequence

def generate_positive_sequences(sequence_length, motifs, motif_frequencies, positional_distribution):
    positive_sequences = []

    for i, motif in enumerate(motifs):
        for seq_count in range(motif_frequencies[i]):
            new_sequence = generate_random_positive_sequence(sequence_length, motif, positional_distribution)
            while new_sequence in positive_sequences:
                new_sequence = generate_random_positive_sequence(sequence_length, motif, positional_distribution)

            positive_sequences.append(new_sequence)

    return positive_sequences

def generate_negative_sequences(sequence_length, motifs, n_negative_sequences, positional_distribution):
    negative_sequences = []

    for i in range(n_negative_sequences):
        new_sequence = generate_random_negative_sequence(sequence_length, motifs, positional_distribution)
        while new_sequence in negative_sequences:
            new_sequence = generate_random_negative_sequence(sequence_length, motifs, positional_distribution)

        negative_sequences.append(new_sequence)

    return negative_sequences


def generate_sequences_with_signal(n_negative_sequences, positional_distribution, motif_frequencies, motif_size, sequence_length):
    motifs = generate_random_motifs(len(motif_frequencies), motif_size, sequence_length, positional_distribution)

    positive_sequences = generate_positive_sequences(sequence_length, motifs, motif_frequencies, positional_distribution)
    negative_sequences = generate_negative_sequences(sequence_length, motifs, n_negative_sequences, positional_distribution)

    return positive_sequences, negative_sequences, motifs

def generate_sequences_without_signal(n_positive_sequences, n_negative_sequences, positional_distribution, sequence_length):

    positive_sequences = generate_negative_sequences(sequence_length, [], n_positive_sequences, positional_distribution)
    negative_sequences = generate_negative_sequences(sequence_length, [], n_negative_sequences, positional_distribution)

    return positive_sequences, negative_sequences, []

def generate_sequences(n_positive_sequences, n_negative_sequences, positional_distribution, motif_frequencies, motif_size, sequence_length):
    if len(motif_frequencies) > 0:
        return generate_sequences_with_signal(n_negative_sequences, positional_distribution, motif_frequencies, motif_size, sequence_length)
    else:
        return generate_sequences_without_signal(n_positive_sequences, n_negative_sequences, positional_distribution, sequence_length)


def export_sequences(file, positive_sequences, negative_sequences):
    df = pd.DataFrame({"file_column_amino_acids": positive_sequences + negative_sequences,
                       "is_binding": [1] * len(positive_sequences) + [0] * len(negative_sequences)})

    df.to_csv(file, index=False)

def compute_motif_frequencies(relative_motif_frequencies, n_positive_sequences):
    return [round(number / sum(relative_motif_frequencies) * n_positive_sequences) for
            number in relative_motif_frequencies]

def parse_positional_distribution_file(positional_distribution_file, sequence_length):
    df = pd.read_csv(positional_distribution_file, sep="\t")

    assert set(df["position"].astype(int)) == set(range(sequence_length)), f"Positions in positional distribution file do not match the given sequence length. Expected positions: {', '.join([str(i) for i in range(sequence_length)])}"

    positional_distribution = {}

    for position in range(sequence_length):
        positional_df = df[df["position"] == position]

        positional_distribution[position] = {"amino_acids": list(positional_df["amino_acid"]),
                                             "relative_frequencies": list(positional_df["relative_frequency"])}

    return positional_distribution

def parse_commandline_arguments(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("--random_seed", default=2022, type=int)
    parser.add_argument("--relative_motif_frequencies", default=[1, 2, 4, 8, 16, 32, 64, 128, 256, 512], type=int, nargs='+')
    parser.add_argument("--motif_size", default=3, type=int)
    parser.add_argument("--n_positive_sequences", default=10000, type=int)
    parser.add_argument("--n_negative_sequences", default=25000, type=int)
    parser.add_argument("--sequence_length", default=10, type=int)

    parser.add_argument("--positional_distribution_file", default="../../data/simulated_data/positional_distribution_mason.tsv", type=str, nargs='+')
    parser.add_argument("--sequences_output_file", default="../../data/simulated_data/sequences.csv", type=str)
    parser.add_argument("--motif_output_file", default="../../data/simulated_data/motifs.tsv", type=str)
    parser.add_argument("--log_file", default="../../data/simulated_data/log.tsv", type=str)

    parsed_args = parser.parse_args(args)
    return parsed_args


if __name__ == "__main__":
    args = parse_commandline_arguments(sys.argv[1:])

    args.motif_frequencies = compute_motif_frequencies(args.relative_motif_frequencies, args.n_positive_sequences)


    build_path(Path(args.log_file).parent)
    build_path(Path(args.sequences_output_file).parent)
    build_path(Path(args.motif_output_file).parent)

    log_arguments(args.log_file, args)

    if args.positional_distribution_file == "":
        positional_distribution = None
    else:
        positional_distribution = parse_positional_distribution_file(args.positional_distribution_file, args.sequence_length)

    random.seed(args.random_seed)
    positive_sequences, negative_sequences, motifs = generate_sequences(n_positive_sequences=args.n_positive_sequences,
                                                                        n_negative_sequences=args.n_negative_sequences,
                                                                        positional_distribution=positional_distribution,
                                                                        motif_frequencies=args.motif_frequencies,
                                                                        motif_size=args.motif_size,
                                                                        sequence_length=args.sequence_length)

    export_sequences(args.sequences_output_file, positive_sequences, negative_sequences)
    write_motifs_to_file(motifs, args.motif_output_file)