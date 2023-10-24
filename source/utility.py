import errno
import os
from pathlib import Path

import numpy as np
import pandas as pd


ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
SEQUENCE_LENGTH = 10

def log_arguments(log_file, args, mode="w"):
    args_str = str(args)[:-1]
    args_str = args_str.replace("Namespace(", "").replace(", ", "\n") + "\n\n"

    with open(log_file, mode) as file:
        file.write(args_str)

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

def get_numpy_sequence_representation(sequence_column, sequence_length=SEQUENCE_LENGTH):
    '''
    Transforms a sequence column (pandas Series) into an optimized numpy representation (2D character array)
    '''
    unicode = sequence_column.to_numpy(dtype=f"U{sequence_length}")
    return unicode.view('U1').reshape(len(sequence_column),-1)


def get_sequences(data, sequence_length, use_np_repr):
    if use_np_repr:
        return get_numpy_sequence_representation(data["file_column_amino_acids"], sequence_length)
    else:
        return data["file_column_amino_acids"]

def read_data_file(data_file, sequence_length=SEQUENCE_LENGTH, use_np_repr=True):
    data = pd.read_csv(data_file, sep=',')

    y_true = data["is_binding"].to_numpy(np.bool_)

    sequences = get_sequences(data, sequence_length, use_np_repr)

    return sequences, y_true


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
