# standard libraries
import unittest

# test library
import ccbbucsd.utilities.analysis_run_prefixes as ns_test


class TestFunctions(unittest.TestCase):
    # No tests for get_* methods, check_or_set, or generate_run_dir because they're so simple

    # region generate_run_prefix
    def test_generate_run_prefix_no_args(self):
        real_output = ns_test.generate_run_prefix(test_val="201700000000")
        self.assertEqual("201700000000", real_output)

    def test_generate_run_prefix_args(self):
        real_output = ns_test.generate_run_prefix("PGP1_MV4", "19mer_1mm", test_val="201700000000")
        self.assertEqual("PGP1_MV4_19mer_1mm_201700000000", real_output)

    # end region

    # region generate_run_prefix_and_dir_dict
    def test_generate_run_prefix_and_dir_dict(self):
        expected_output = {'run_dir': '/my/testdir/PGP1_MV4_19mer_1mm_201700000000',
                           'run_prefix': 'PGP1_MV4_19mer_1mm_201700000000'}

        real_output = ns_test.generate_run_prefix_and_dir_dict("/my/testdir/", "PGP1_MV4", "19mer_1mm",
                                                               test_val="201700000000")
        self.assertEqual(expected_output, real_output)

    # end region

    # region strip_run_prefix
    def test_strip_run_prefix_w_separator_after(self):
        real_output = ns_test.strip_run_prefix("PGP1_MV4_19mer_1mm_201700000000_counts_combined.txt",
                                               "PGP1_MV4_19mer_1mm_201700000000")
        self.assertEqual("counts_combined.txt", real_output)

    def test_strip_run_prefix_w_separator_before(self):
        real_output = ns_test.strip_run_prefix("counts_combined_PGP1_MV4_19mer_1mm_201700000000.txt",
                                               "PGP1_MV4_19mer_1mm_201700000000")
        self.assertEqual("counts_combined.txt", real_output)

    def test_strip_run_prefix_without_separator(self):
        real_output = ns_test.strip_run_prefix("counts_combinedPGP1_MV4_19mer_1mm_201700000000.txt",
                                               "PGP1_MV4_19mer_1mm_201700000000")
        self.assertEqual("counts_combined.txt", real_output)

    def test_strip_run_prefix_error(self):
        with self.assertRaises(ValueError):
            real_output = ns_test.strip_run_prefix("    _PGP1_MV4_19mer_1mm_201700000000 ",
                                                   "PGP1_MV4_19mer_1mm_201700000000")

    # end region

    # region describe_var_list
    def test_describe_var_list(self):
        expected_output = '__author__: Amanda Birmingham\n__maintainer__: Amanda Birmingham\n'
        real_output = ns_test.describe_var_list(["__author__", "__maintainer__"])
        self.assertEqual(expected_output, real_output)

    # end region

