import errno
import os
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.metrics import recall_score, f1_score, precision_score, accuracy_score, balanced_accuracy_score

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
SEQUENCE_LENGTH = 10

def log_arguments(log_file, args, mode="w"):
    args_str = str(args)[:-1]
    args_str = args_str.replace("Namespace(", "").replace(", ", "\n") + "\n\n"

    with open(log_file, mode) as file:
        file.write(args_str)

# first version: using pandas
def pd_test_aa(sequence_column, index, aa):
    return (sequence_column.str.get(index) == aa).to_numpy()

def pd_test_position(sequence_column, index, aas):
    return np.logical_or.reduce([pd_test_aa(sequence_column, index, aa) for aa in aas])

# optimized version: using numpy
def np_test_aa(sequences, index, aa):
    return sequences[:, index] == aa

def np_test_position(sequences, index, aas):
    return np.logical_or.reduce([np_test_aa(sequences, index, aa) for aa in aas])

def test_motif(sequences, indices, amino_acids, test_fn=np_test_position):
    '''
    Tests for all sequences whether it contains the given motif (defined by indices and amino acids)

    test_fn: should be np_test_aa/pd_test_aa (when testing a single amino acid for each position) or pd_test_position/np_test_position (when multiple amino acids may be OR'd at one position)
    '''
    return np.logical_and.reduce([test_fn(sequences, index, amino_acid) for index, amino_acid in zip(indices, amino_acids)])

def test_rule_tree(sequences, rule_tree):
    return np.logical_or.reduce([test_motif(sequences, indices, amino_acids) for indices, amino_acids in rule_tree])

def performance_rule_tree(sequences, y_true, weights, rule_tree, performance_metric=accuracy_score):
    pred = test_rule_tree(sequences, rule_tree)
    return performance_metric(y_true=y_true, y_pred=pred, sample_weight=weights)

def compute_column_contributions(column, alphabet=ALPHABET, pseudocount_value=1):
    # if weighting is to be included here, then that weight should be applied to aa_list.count(amino_acid)
    aa_list = list(column)
    total = len(aa_list) + pseudocount_value

    normalized_count = lambda amino_acid: (aa_list.count(amino_acid) + pseudocount_value) / total
    return {amino_acid: normalized_count(amino_acid) for amino_acid in alphabet}

def compute_column_frequencies(column, pseudocount_value=1, alphabet=ALPHABET):
    aa_list = list(column)
    return {amino_acid: aa_list.count(amino_acid) + pseudocount_value for amino_acid in alphabet}

def compute_positional_aa_contributions(sequence_column, sequence_length=SEQUENCE_LENGTH, pseudocount_value=1):
    return {position: compute_column_contributions(sequence_column.str.get(position), pseudocount_value=pseudocount_value) for position in range(sequence_length)}

def compute_raw_sequence_weight(sequence, positional_weights, alphabet=ALPHABET):
    sequence_weight = 1

    uniform = 1 / len(alphabet)

    for i, aa in enumerate(sequence):
        sequence_weight *= uniform / positional_weights[i][aa]

    return sequence_weight

def compute_sequence_weights_train_test_split(train_sequence_column, test_sequence_column, sequence_length, weight_threshold, alphabet=ALPHABET, pseudocount_value=1):
    positional_weights = compute_positional_aa_contributions(train_sequence_column, sequence_length, pseudocount_value=pseudocount_value)
    weights = np.array([compute_raw_sequence_weight(sequence, positional_weights, alphabet=alphabet) for sequence in test_sequence_column])
    weights[weights > weight_threshold] = weight_threshold
    return weights

def compute_sequence_weights(sequence_column, sequence_length, weight_threshold, alphabet=ALPHABET, pseudocount_value=1):
    return compute_sequence_weights_train_test_split(sequence_column, sequence_column, sequence_length, weight_threshold, alphabet, pseudocount_value)

def get_numpy_sequence_representation(sequence_column, sequence_length=SEQUENCE_LENGTH):
    '''
    Transforms a sequence column (pandas Series) into an optimized numpy representation (2D character array)
    '''
    unicode = sequence_column.to_numpy(dtype=f"U{sequence_length}")
    return unicode.view('U1').reshape(len(sequence_column),-1)

def count_true_positives(y_true, y_pred):
    return sum(np.logical_and(y_true, y_pred))

def get_indices_from_str(indices_str):
    return tuple([int(i) for i in indices_str])

def get_sequences(data, sequence_length, use_np_repr):
    if use_np_repr:
        return get_numpy_sequence_representation(data["file_column_amino_acids"], sequence_length)
    else:
        return data["file_column_amino_acids"]

def read_data_file(data_file, weight_threshold=np.inf, sequence_length=SEQUENCE_LENGTH, alphabet=ALPHABET, use_np_repr=True,
                   calculate_weights = True):
    data = pd.read_csv(data_file, sep=',')

    if calculate_weights:
        weights = compute_sequence_weights(data["file_column_amino_acids"], sequence_length, weight_threshold, alphabet=alphabet)
    else:
        weights = None

    y_true = data["is_binding"].to_numpy(np.bool_)

    sequences = get_sequences(data, sequence_length, use_np_repr)

    return sequences, y_true, weights

def read_test_data_file(train_data_file, test_data_file, weight_threshold, sequence_length=SEQUENCE_LENGTH, alphabet=ALPHABET, use_np_repr=True, calculate_weights = True):
    # todo generalize with compute_sequence_weights func! give in positional_weights as arg!
    train_data = pd.read_csv(train_data_file, sep=',')
    test_data = pd.read_csv(test_data_file, sep=',')

    if calculate_weights:
        # compute weights using train data positional weights
        weights = compute_sequence_weights_train_test_split(train_data["file_column_amino_acids"], test_data["file_column_amino_acids"], sequence_length,
                                                            weight_threshold, alphabet=alphabet)
    else:
        weights = None

    y_true = test_data["is_binding"].to_numpy(np.bool_)

    sequences = get_sequences(test_data, sequence_length, use_np_repr)

    return sequences, y_true, weights


def read_positions_file(positions_file):
    df = pd.read_csv(positions_file, sep="\t", dtype=str)

    df["indices"] = df["indices"].apply(get_indices_from_str)
    df["amino_acids"] = df["amino_acids"].str.split(",")
    df["positive_predictions"] = df["positive_predictions"].astype(int)
    df["raw_tp"] = df["raw_tp"].astype(int)

    for colname in ["raw_precision", "w_precision", "raw_recall", "w_recall", "raw_f1", "w_f1"]:
        df[colname] = df[colname].astype(float)

    return df

def read_positions_files(positions_files):
    return pd.concat([read_positions_file(file) for file in positions_files], ignore_index=True)

def report_stats(y_true, pred, weights):
    print("weighted:")
    print("- f1", f1_score(y_true=y_true, y_pred=pred, sample_weight=weights))
    print("- precision", precision_score(y_true=y_true, y_pred=pred, sample_weight=weights))
    print("- recall", recall_score(y_true=y_true, y_pred=pred, sample_weight=weights))
    print("- accuracy", accuracy_score(y_true=y_true, y_pred=pred, sample_weight=weights))
    print("- balanced accuracy", balanced_accuracy_score(y_true=y_true, y_pred=pred, sample_weight=weights))

    print("\nunweighted:")
    print("- f1", f1_score(y_true=y_true, y_pred=pred))
    print("- precision", precision_score(y_true=y_true, y_pred=pred))
    print("- recall", recall_score(y_true=y_true, y_pred=pred))
    print("- accuracy", accuracy_score(y_true=y_true, y_pred=pred))
    print("- balanced accuracy", balanced_accuracy_score(y_true=y_true, y_pred=pred))

def report_stats_wide(y_true, pred, weights):
    print("w_f1\tw_precision\tw_recall\tw_accuracy\tw_balanced_accuracy\tf1\tprecision\trecall\taccuracy\tbalanced_accuracy")
    print("\t".join([str(score) for score in
                     [f1_score(y_true=y_true, y_pred=pred, sample_weight=weights),
                     precision_score(y_true=y_true, y_pred=pred, sample_weight=weights),
                     recall_score(y_true=y_true, y_pred=pred, sample_weight=weights),
                     accuracy_score(y_true=y_true, y_pred=pred, sample_weight=weights),
                     balanced_accuracy_score(y_true=y_true, y_pred=pred, sample_weight=weights),
                     f1_score(y_true=y_true, y_pred=pred),
                     precision_score(y_true=y_true, y_pred=pred),
                     recall_score(y_true=y_true, y_pred=pred),
                     accuracy_score(y_true=y_true, y_pred=pred),
                     balanced_accuracy_score(y_true=y_true, y_pred=pred)]]))


def print_rule_tree(rule_tree):
    for rule in rule_tree:
        print("(" + " ".join([str(i) for i in rule[0]]) + ") (" +  " ".join(rule[1]) + ")")


def motif_to_string(indices, amino_acids, value_sep="&", motif_sep="\t", newline=True):
    suffix = "\n" if newline else ""
    return f"{value_sep.join([str(idx) for idx in indices])}{motif_sep}{value_sep.join(amino_acids)}{suffix}"


def string_to_motif(string, value_sep, motif_sep):
    indices_str, amino_acids_str = string.strip().split(motif_sep)
    indices = [int(i) for i in indices_str.split(value_sep)]
    amino_acids = amino_acids_str.split(value_sep)
    return indices, amino_acids


def read_motifs_from_file(filepath):
    with open(filepath) as file:
        header = file.readline()
        assert header == "indices\tamino_acids\n", f"Motif file at {filepath} is expected to contain this header: " \
                                                   f"'indices\tamino_acids', found the following instead: '{header}'"

        motifs = [string_to_motif(line, value_sep="&", motif_sep="\t") for line in
                  file.readlines()]

    return motifs

def filter_single_motifs(motifs):
    return list(filter(lambda motif: all([len(aa) == 1 for aa in motif[1]]), motifs))

def get_confusion_matrix_subset(np_sequences, y_true, motif, group):
    contains_motif = test_motif(np_sequences, motif[0], motif[1])

    if group == "P":
        is_motif_sequence = contains_motif
    elif group == "TP":
        is_motif_sequence = contains_motif & y_true
    elif group == "FP":
        is_motif_sequence = contains_motif & ~y_true
    elif group == "N":
        is_motif_sequence = ~contains_motif
    elif group == "TN":
        is_motif_sequence = ~contains_motif & ~y_true
    elif group == "FN":
        is_motif_sequence = ~contains_motif & y_true
    else:
        raise ValueError(f"Illegal motif_sequence_group: {group}")

    return np_sequences[is_motif_sequence]

def get_color_discrete_sequence():
    import plotly.express as px
    return px.colors.qualitative.Pastel[:-1] + px.colors.qualitative.Set3

def motif_to_string(indices, amino_acids, value_sep="&", motif_sep="\t", newline=True):
    suffix = "\n" if newline else ""
    return f"{value_sep.join([str(idx) for idx in indices])}{motif_sep}{value_sep.join(amino_acids)}{suffix}"

def write_motifs_to_file(motifs, filepath):
    with open(filepath, "w") as file:
        file.write("indices\tamino_acids\n")

        for indices, amino_acids in motifs:
            file.write(motif_to_string(indices, amino_acids))

def build_path(path):
    path = Path(path)
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    return path
