'''
Code to create Supplementary Table 6:
The number of motif-containing sequences in the simulated dataset (training, validation and test subsets).
'''

from source.utility import *


motifs = [([4, 7, 8],["G", "P", "I"]),
          ([0, 5, 6],["Y", "M", "F"]),
          ([0, 7, 8],["Y", "V", "F"]),
          ([4, 5, 6],["S", "S", "F"]),
          ([0, 2, 8],["F", "L", "F"]),
          ([1, 6, 9],["N", "Y", "V"]),
          ([5, 6, 8],["S", "Y", "F"]),
          ([0, 2, 8],["W", "P", "L"]),
          ([5, 6, 8],["F", "Y", "L"]),
          ([1, 6, 8],["R", "Y", "N"])]


train, true_train = read_data_file(f"../../data/simulated_data/train_val_test/training.csv")
val, true_val = read_data_file(f"../../data/simulated_data/train_val_test/validation.csv")
test, true_test = read_data_file(f"../../data/simulated_data/train_val_test/test.csv")
all, true_all = read_data_file(f"../../data/simulated_data/train_val_test/training_validation_test.csv")

train = train[true_train]
val = val[true_val]
test = test[true_test]
all = all[true_all]

def reformat_motif(motif, motif_size=3):
    return ", ".join([f"{motif[1][i]}{motif[0][i]}" for i in range(motif_size)])


print("motif\ttotal\ttrain\tval\ttest")

for i, motif in enumerate(motifs):

    n_total = sum(test_motif(all, motif[0], motif[1]))
    n_train = sum(test_motif(train, motif[0], motif[1]))
    n_val = sum(test_motif(val, motif[0], motif[1]))
    n_test= sum(test_motif(test, motif[0], motif[1]))

    print(f"{reformat_motif(motif)}\t{n_total}\t{n_train}\t{n_val}\t{n_test}")
