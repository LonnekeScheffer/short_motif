import pandas as pd
from unittest import TestCase
from source.utility import *

class TestIdentifyHighPrecisionMotifs(TestCase):
    def test_pd_test_aa(self):
        sequence_column = pd.Series(["AA", "BB", "AB", "BA"])

        self.assertListEqual(list(pd_test_aa(sequence_column, 0, "A")), [True, False, True, False])
        self.assertListEqual(list(pd_test_aa(sequence_column, 1, "A")), [True, False, False, True])
        self.assertListEqual(list(pd_test_aa(sequence_column, 0, "B")), [False, True, False, True])
        self.assertListEqual(list(pd_test_aa(sequence_column, 1, "B")), [False, True, True, False])

    def test_pd_test_aa(self):
        sequence_array = np.asarray(['A' 'A', 'B' 'B', 'A' 'B', 'B' 'A']).view('U1').reshape(4, -1)

        self.assertListEqual(list(np_test_aa(sequence_array, 0, "A")), [True, False, True, False])
        self.assertListEqual(list(np_test_aa(sequence_array, 1, "A")), [True, False, False, True])
        self.assertListEqual(list(np_test_aa(sequence_array, 0, "B")), [False, True, False, True])
        self.assertListEqual(list(np_test_aa(sequence_array, 1, "B")), [False, True, True, False])

    def test_test_position(self):
        sequence_array = np.asarray(['A' 'A', 'B' 'B', 'A' 'B', 'C' 'A']).view('U1').reshape(4, -1)

        self.assertListEqual(list(np_test_position(sequence_array, 0, "A")), [True, False, True, False])
        self.assertListEqual(list(np_test_position(sequence_array, 0, "AB")), [True, True, True, False])
        self.assertListEqual(list(np_test_position(sequence_array, 0, "BC")), [False, True, False, True])
        self.assertListEqual(list(np_test_position(sequence_array, 0, "ABC")), [True, True, True, True])

    def test_test_motif(self):
        sequence_array = np.asarray(['A' 'A', 'B' 'B', 'A' 'B', 'B' 'A']).view('U1').reshape(4, -1)

        self.assertListEqual(list(test_motif(sequence_array, (0, 1), ("A", "B"))), [False, False, True, False])
        self.assertListEqual(list(test_motif(sequence_array, (0, 1), ("E", "E"))), [False, False, False, False])
        self.assertListEqual(list(test_motif(sequence_array, (0, 1), ("DE", "DE"))), [False, False, False, False])
        self.assertListEqual(list(test_motif(sequence_array, (0, 1), ("A", "BA"))), [True, False, True, False])
        self.assertListEqual(list(test_motif(sequence_array, (0, 1), ("AB", "AB"))), [True, True, True, True])
        self.assertListEqual(list(test_motif(sequence_array, (0, 1), ("C", "AB"))), [False, False, False, False])
        self.assertListEqual(list(test_motif(sequence_array, (0, 1), ("AB", "C"))), [False, False, False, False])

    def test_get_numpy_sequence_representation(self):
        sequence_column = pd.Series(["AA", "BB", "AB", "BA"])
        output = get_numpy_sequence_representation(sequence_column, 2)

        expected = np.asarray(['A' 'A', 'B' 'B', 'A' 'B', 'B' 'A']).view('U1').reshape(4, -1)

        self.assertEqual(output.shape, expected.shape)

        for i in range(len(output)):
            self.assertListEqual(list(output[i]), list(expected[i]))

            for j in range(len(output[i])):
                self.assertEqual(type(output[i][j]), type(expected[i][j]))
