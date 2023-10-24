import itertools
from math import prod

import pandas as pd

distribution = pd.read_csv("../../data/simulated_data/positional_distribution_mason.tsv", sep="\t")
distribution = distribution[distribution["relative_frequency"] > 0]

print(distribution)
n_options_per_position = list(distribution.position.value_counts())

print(f"Number of different possible amino acids per sequence position: {n_options_per_position}")


for motif_size in range(1, 6):
    n_motifs = 0

    print(f"Motif size: {motif_size}")
    for position_combo in itertools.combinations(n_options_per_position, r=motif_size):
        n_motifs += prod(position_combo)

    print(f"Theoretically possible motifs: {n_motifs}")
