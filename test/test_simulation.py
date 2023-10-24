from unittest import TestCase
from source.simulation.simulate_synthetic_data import *

class TestSimulation(TestCase):

    def test_get_amino_acids_from_positional_distribution(self):
        result = get_amino_acids_from_positional_distribution([1, 3], {1: {"amino_acids": ["a", "b"],
                                                                           "relative_frequencies": [1, 0]},
                                                                       2: {"amino_acids": ["c", "d"],
                                                                           "relative_frequencies": [1, 0]},
                                                                       3: {"amino_acids": ["e", "f"],
                                                                           "relative_frequencies": [0, 1]}})

        self.assertEqual(result, ["a", "f"])

    def test_generate_random_motif(self):
        result = generate_random_motif(1, sequence_length=3, positional_distribution={0: {"amino_acids": ["a", "b"],
                                                                                          "relative_frequencies": [1, 0]},
                                                                                      1: {"amino_acids": ["a", "b"],
                                                                                          "relative_frequencies": [1, 0]},
                                                                                      2: {"amino_acids": ["a", "b"],
                                                                                          "relative_frequencies": [1, 0]}})

        self.assertEqual(result[1], ("a",))
        self.assertTrue(result[0]==(0,) or result[0]==(1,) or result[0]==(2,))

    def test_generate_random_positive_sequence(self):
        result = generate_random_positive_sequence(sequence_length=3,
                                                   implant_motif=[(0,), ("b",)],
                                                   positional_distribution={0: {"amino_acids": ["a", "b"],
                                                                                "relative_frequencies": [1, 0]},
                                                                            1: {"amino_acids": ["a", "b"],
                                                                                "relative_frequencies": [1, 0]},
                                                                            2: {"amino_acids": ["a", "b"],
                                                                                "relative_frequencies": [1, 0]}})

        self.assertEqual(result, "baa")

    def test_test_contains_any_motif(self):
        self.assertTrue(test_contains_any_motif("aaa", [((0,), ("a",)), ((0,), ("b",))]))
        self.assertTrue(test_contains_any_motif("aaa", [[(0,), ("a",)]]))
        self.assertFalse(test_contains_any_motif("aaa", [[(0,), ("b",)]]))

    def test_generate_random_negative_sequence(self):
        result = generate_random_negative_sequence(3, [[(0,), ("a",)]],
                                                   positional_distribution={0: {"amino_acids": ["a", "b"],
                                                                                "relative_frequencies": [10, 1]},
                                                                            1: {"amino_acids": ["a", "b"],
                                                                                "relative_frequencies": [1, 0]},
                                                                            2: {"amino_acids": ["a", "b"],
                                                                                "relative_frequencies": [1, 0]}})

        self.assertEqual(result, "baa")

    def test_generate_positive_sequences(self):
        result = generate_positive_sequences(sequence_length=2,
                                             motifs=[[(0,), ("a",)]],
                                             motif_frequencies=[5],
                                             positional_distribution={0: {"amino_acids": ["a", "b"],
                                                                          "relative_frequencies": [1, 1]},
                                                                      1: {"amino_acids": ["a", "b", "c", "d", "e"],
                                                                          "relative_frequencies": [1, 1, 1, 1, 1]}})

        self.assertListEqual(sorted(result), sorted(["aa", "ab", "ac", "ad", "ae"]))

    def test_generate_negative_sequences(self):
        result = generate_negative_sequences(sequence_length=2,
                                             motifs=[[(0,), ("a",)]],
                                             n_negative_sequences=5,
                                             positional_distribution={0: {"amino_acids": ["a", "b"],
                                                                          "relative_frequencies": [1, 1]},
                                                                      1: {"amino_acids": ["a", "b", "c", "d", "e"],
                                                                          "relative_frequencies": [1, 1, 1, 1, 1]}})

        self.assertListEqual(sorted(result), sorted(["ba", "bb", "bc", "bd", "be"]))

    def test_compute_motif_frequencies(self):
        result = compute_motif_frequencies([1, 1, 2], 20)

        self.assertEqual(result, [5, 5, 10])
