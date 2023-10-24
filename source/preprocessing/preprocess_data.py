import errno

import os
from source.preprocessing.train_val_test_split import *



def build_path(path):
    path = Path(path)
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    return path

def aggregate_duplicates_within_class(data):
    '''
    Remove duplicate amino acid sequences (different nucleotide sequences),
    sum the counts of duplicates
    '''

    data["count"] = data.groupby(["file_column_amino_acids"])["count"].transform("sum")
    data = data[["count", "file_column_amino_acids", "is_binding"]].copy()
    data.drop_duplicates(inplace=True, keep="first")

    return data

def annotate_with_ratios(duplicates):
    binder_duplicates = duplicates[duplicates["is_binding"] == 1][["count", "file_column_amino_acids"]].copy()
    nonbinder_duplicates = duplicates[duplicates["is_binding"] == 0][["count", "file_column_amino_acids"]].copy()

    binder_duplicates.rename(columns={"count": "count_1"}, inplace=True)
    nonbinder_duplicates.rename(columns={"count": "count_0"}, inplace=True)

    duplicates = pd.merge(binder_duplicates, nonbinder_duplicates, on="file_column_amino_acids")
    duplicates["01_ratio"] = duplicates["count_0"] / duplicates["count_1"]
    duplicates["10_ratio"] = duplicates["count_1"] / duplicates["count_0"]

    return duplicates

def select_best_class(duplicates, min_ratio):
    '''only keep sequences where one of the classes has at least min_ratio x as high count as the other class'''
    duplicates = duplicates[duplicates[["01_ratio", "10_ratio"]].max(axis=1) > min_ratio].copy()
    duplicates["count"] = duplicates[["count_0", "count_1"]].max(axis=1)
    duplicates["alt_count"] = duplicates[["count_0", "count_1"]].min(axis=1)

    duplicates["is_binding"] = (duplicates["count_1"] > duplicates["count_0"]).astype(int)

    return duplicates

def select_duplicates_across_classes(data, min_ratio=2):
    '''
    For sequences appearing in both binder and nonbinder classes, select the class
    with the highest count (only if highest count is greater than lowest count * min_ratio)
    '''
    duplicate_indices = data.duplicated(subset=["file_column_amino_acids"], keep=False)

    unique = data[~duplicate_indices].copy()
    unique["alt_count"] = 0

    duplicates = data[duplicate_indices].copy()
    duplicates = annotate_with_ratios(duplicates)
    duplicates = select_best_class(duplicates, min_ratio)

    duplicates = duplicates[unique.columns].copy()

    if len(duplicates) > 0:
        data = pd.concat([unique, duplicates])
    else:
        data = unique

    data.sort_values(by="is_binding", inplace=True)
    return data

def combine_dfs(positives_df, negatives_df, min_ratio=2):
    positives_df = aggregate_duplicates_within_class(positives_df)
    negatives_df = aggregate_duplicates_within_class(negatives_df)

    data = pd.concat([positives_df, negatives_df])
    data = select_duplicates_across_classes(data, min_ratio)

    return data.reindex(columns=["file_column_amino_acids", "is_binding", "count", "alt_count"])


def export_files(data, output_folder, dataset_name):
    build_path(output_folder)
    full_dataset_path = str(Path(output_folder) / f"{dataset_name}_all.csv")

    data.to_csv(full_dataset_path, index=False)
    args = parse_data_splitting_arguments(["--input_file", full_dataset_path, "--output_folder", output_folder])
    write_data_splits(args)


def get_mason_df(positives_file_path, negatives_file_path):
    column_mapping = {"Count": "count", "AASeq": "file_column_amino_acids", "AgClass": "is_binding"}

    positives_df = pd.read_csv(positives_file_path).rename(columns=column_mapping, inplace=False)
    negatives_df = pd.read_csv(negatives_file_path).rename(columns=column_mapping, inplace=False)

    return combine_dfs(positives_df, negatives_df)


def get_mehta_df(positives_file_path, negatives_file_path):
    column_mapping = {"aa_cdr3_seq": "file_column_amino_acids"}

    positives_df = pd.read_csv(positives_file_path).rename(columns=column_mapping, inplace=False)
    negatives_df = pd.read_csv(negatives_file_path).rename(columns=column_mapping, inplace=False)

    positives_df["is_binding"] = 1
    negatives_df["is_binding"] = 0

    return combine_dfs(positives_df, negatives_df)


if __name__ == "__main__":
    data = get_mason_df(positives_file_path="../../data/mason/from_github/mHER_H3_AgPos.csv",
                        negatives_file_path="../../data/mason/from_github/mHER_H3_AgNeg.csv")

    export_files(data,
                 output_folder=f"../../data/mason/train_val_test",
                 dataset_name="mason")


