# standard libraries
import unittest

# ccbb libraries
from ccbbucsd.utilities.bio_seq_utilities import rev_comp_canonical_dna_seq, trim_seq

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
    # region trim_seq tests
    def test_trim_seq_long(self):
        input_seq = "ACGT"
        retain_len = 3

        # trim from 5p end
        output_5p = trim_seq(input_seq, retain_len, False)
        self.assertEqual("CGT", output_5p)

        # trim from 3p end
        output_3p = trim_seq(input_seq, retain_len, True)
        self.assertEqual("ACG", output_3p)

    def test_trim_seq_exact(self):
        input_seq = "ACGT"
        retain_len = 4

        # trim from 5p end
        output_5p = trim_seq(input_seq, retain_len, False)
        self.assertEqual(input_seq, output_5p)

        # trim from 3p end
        output_3p = trim_seq(input_seq, retain_len, True)
        self.assertEqual(input_seq, output_3p)

    def test_trim_seq_short(self):
        input_seq = "ACGT"
        retain_len = 5

        # trim from 5p end
        with self.assertRaises(ValueError):
            trim_seq(input_seq, retain_len, False)

        # trim from 3p end
        with self.assertRaises(ValueError):
            trim_seq(input_seq, retain_len, True)

    # endregion

    # region rev_comp_canonical_dna_seq tests
    def test_rev_comp_canonical_dna_seq(self):
        input = "ACGccAgtaT"
        expected_output = "AtacTggCGT"

        real_output = rev_comp_canonical_dna_seq(input)
        self.assertEqual(expected_output, real_output)
    # endregion
