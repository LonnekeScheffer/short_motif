import sys
import argparse
from pathlib import Path
import pandas as pd
import random
from source.utility import build_path

def add_identifiers_column(df):
    df["sequence_identifiers"] = list(range(len(df)))

def get_train_val_test_indices(n_examples, args):
    assert args.training_percentage + args.val_percentage + args.test_percentage == 1, f"{args.training_percentage + args.val_percentage + args.test_percentage} != 1"

    indices = list(range(n_examples))

    random.seed(args.random_seed)
    random.shuffle(indices)

    limit = int(n_examples * args.training_percentage)
    train_indices = indices[:limit]
    val_test_indices = indices[limit:]

    limit = int(n_examples * args.val_percentage)
    val_indices = val_test_indices[:limit]
    test_indices = val_test_indices[limit:]

    return sorted(train_indices), sorted(val_indices), sorted(test_indices)


def write_data_splits(args):
    build_path(args.output_folder)

    df = pd.read_csv(args.input_file)

    if args.add_identifiers:
        add_identifiers_column(df)

    train_indices, val_indices, test_indices = get_train_val_test_indices(len(df), args)

    training_df = df.iloc[train_indices].copy()
    val_df = df.iloc[val_indices].copy()
    test_df = df.iloc[test_indices].copy()
    train_val_df = df.iloc[train_indices + val_indices].copy()

    df.to_csv(Path(args.output_folder) / "training_validation_test.csv", sep=",", index=False)
    training_df.to_csv(Path(args.output_folder) / "training.csv", sep=",", index=False)
    val_df.to_csv(Path(args.output_folder) / "validation.csv", sep=",", index=False)
    test_df.to_csv(Path(args.output_folder) / "test.csv", sep=",", index=False)
    train_val_df.to_csv(Path(args.output_folder) / "training_validation.csv", sep=",", index=False)

    df.rename(columns={"sequence_identifiers": "example_id"}, inplace=True)
    training_df.rename(columns={"sequence_identifiers": "example_id"}, inplace=True)
    val_df.rename(columns={"sequence_identifiers": "example_id"}, inplace=True)
    test_df.rename(columns={"sequence_identifiers": "example_id"}, inplace=True)
    train_val_df.rename(columns={"sequence_identifiers": "example_id"}, inplace=True)

    if args.add_identifiers:
        df["example_id"].to_csv(Path(args.output_folder) / "training_validation_test_identifiers.txt", index=False)
        training_df["example_id"].to_csv(Path(args.output_folder) / "training_identifiers.txt", index=False)
        val_df["example_id"].to_csv(Path(args.output_folder) / "validation_identifiers.txt", index=False)
        test_df["example_id"].to_csv(Path(args.output_folder) / "test_identifiers.txt", index=False)
        train_val_df["example_id"].to_csv(Path(args.output_folder) / "training_validation_identifiers.txt", index=False)


def parse_data_splitting_arguments(args):
    parser = argparse.ArgumentParser()

    # parser.add_argument("--input_file", default="/Users/lonneke/PycharmProjects/rule_based_classifier/data/mason/train_val_test/mason_all.csv")
    # parser.add_argument("--output_folder", default="/Users/lonneke/PycharmProjects/rule_based_classifier/data/mason/train_val_test")
    # parser.add_argument("--input_file", default="/Users/lonneke/PycharmProjects/rule_based_classifier/data/simulated_data/skewed_simulation/sequences.csv")
    # parser.add_argument("--output_folder", default="/Users/lonneke/PycharmProjects/rule_based_classifier/data/simulated_data/skewed_simulation/train_val_test")
    parser.add_argument("--input_file", default="/Users/lonneke/PycharmProjects/rule_based_classifier/data/simulated_data/no_motifs/sequences.csv")
    parser.add_argument("--output_folder", default="/Users/lonneke/PycharmProjects/rule_based_classifier/data/simulated_data/no_motifs/train_val_test")

    parser.add_argument("--add_identifiers", default=True, type=bool)
    parser.add_argument("--training_percentage", default=0.5, type=float)
    parser.add_argument("--val_percentage", default=0.25, type=float)
    parser.add_argument("--test_percentage", default=0.25, type=float)

    parser.add_argument("--random_seed", default=2022, type=int)

    parsed_args = parser.parse_args(args)

    return parsed_args


if __name__ == "__main__":
    args = parse_data_splitting_arguments(sys.argv[1:])
    write_data_splits(args)
