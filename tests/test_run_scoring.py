# standard libraries
import tempfile
import unittest
import warnings

from dual_crispr import run_scoring as ns_test


class TestFunctions(unittest.TestCase):
    # no tests for _parse_cmd_line_args as it is so simple
    # no tests for main as it just chains together calls to other tested methods

    # region _set_params
    def test__set_params_test(self):
        config_str = """[DEFAULT]
machine_configuration = laptop

[c4_2xlarge]
num_processors: 7

[laptop]
num_processors: 3

[count_pipeline]
keep_gzs: False
full_5p_r1: TATATATCTTGTGGAAAGGACGAAACACCG
full_5p_r2: CCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC
full_3p_r1: GTTTCAGAGCTATGCTGGAAACTGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAG
full_3p_r2: CAAACAAGGCTTTTCTCCAAGGGATATTTATAGTCTCAAAACACACAATTACTTTACAGTTAGGGTGAGTTTCCTTTTGTGCTGTTTTTTAAAATA
len_of_seq_to_match = 19
num_allowed_mismatches = 1
# Set-up for pipeline notebooks; do not modify unless you are a power user!
notebook_basenames_list: Dual CRISPR 1-Construct Scaffold Trimming.ipynb,Dual CRISPR 2-Constuct Filter.ipynb,Dual CRISPR 3-Construct Counting.ipynb,Dual CRISPR 4-Count Combination.ipynb,Dual CRISPR 5-Count Plots.ipynb

[score_pipeline]
time_prefixes: T,D
# min_count_limit in absolute counts, not log2
min_count_limit: 10
# max_fraction_acceptable_spline_density_diff is % of diff between max spline and min density
max_fraction_acceptable_spline_density_diff: 0.02
# any threshold throwing out > max_fraction_counts_excluded% of counts is not acceptable
max_fraction_counts_excluded: 0.95
use_seed = False
num_iterations = 1000
# Set-up for pipeline notebooks; do not modify unless you are a power user!
notebook_basenames_list: Dual CRISPR 6-Scoring Preparation.ipynb,Dual CRISPR 7-Abundance Thresholds.ipynb,Dual CRISPR 8-Construct Scoring.ipynb

[test]
use_seed = True
num_iterations = 2
"""

        count_fps_or_dirs = "/my/counts/dir1,/my/counts_dir2"
        day_timepoints_str = "3,10,21"
        provided_output_dir = "/my/output_parent"
        is_test = True

        # count_fps_or_dirs string, day_timepoints_str, time_prefixes, and notebooks list should NOT be parsed to list
        # use_seed SHOULD be parsed to boolean
        # num_iterations SHOULD be parsed to int
        expected_output = {'machine_configuration': "laptop",
                           'count_fps_or_dirs': '/my/counts/dir1,/my/counts_dir2',
                           'day_timepoints_str': '3,10,21',
                           'processed_data_dir': "/my/output_parent",
                           'notebook_basenames_list': 'Dual CRISPR 6-Scoring Preparation.ipynb,Dual '
                            'CRISPR 7-Abundance Thresholds.ipynb,Dual CRISPR '
                             '8-Construct Scoring.ipynb',
                           'min_count_limit': 10,
                           'max_fraction_acceptable_spline_density_diff': 0.02,
                           'max_fraction_counts_excluded': 0.95,
                           'num_iterations': 2,
                           'time_prefixes': 'T,D',
                           'use_seed': True}

        temp_config = tempfile.NamedTemporaryFile(mode="w")
        temp_config.write(config_str)
        temp_config.seek(0)

        with warnings.catch_warnings(record=True) as warnings_list:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")

            real_output = ns_test._set_params(count_fps_or_dirs, day_timepoints_str, provided_output_dir, is_test,
                                              temp_config.name)

            assert len(warnings_list) == 1
            assert "Scoring is running in TEST MODE; do not use results for data analysis!" in str(warnings_list[-1].message)

        self.assertEqual(expected_output, real_output)

    def test_set_params_real(self):
        config_str = """[DEFAULT]
machine_configuration = laptop

[c4_2xlarge]
num_processors: 7
keep_gzs: False

[laptop]
num_processors: 3
keep_gzs: True


[count_pipeline]
full_5p_r1: TATATATCTTGTGGAAAGGACGAAACACCG
full_5p_r2: CCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC
full_3p_r1: GTTTCAGAGCTATGCTGGAAACTGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAG
full_3p_r2: CAAACAAGGCTTTTCTCCAAGGGATATTTATAGTCTCAAAACACACAATTACTTTACAGTTAGGGTGAGTTTCCTTTTGTGCTGTTTTTTAAAATA
len_of_seq_to_match = 19
num_allowed_mismatches = 1
# Set-up for pipeline notebooks; do not modify unless you are a power user!
notebook_basenames_list: Dual CRISPR 1-Construct Scaffold Trimming.ipynb,Dual CRISPR 2-Constuct Filter.ipynb,Dual CRISPR 3-Construct Counting.ipynb,Dual CRISPR 4-Count Combination.ipynb,Dual CRISPR 5-Count Plots.ipynb

[score_pipeline]
time_prefixes: T,D
# min_count_limit in absolute counts, not log2
min_count_limit: 10
# max_fraction_acceptable_spline_density_diff is % of diff between max spline and min density
max_fraction_acceptable_spline_density_diff: 0.02
# any threshold throwing out > max_fraction_counts_excluded% of counts is not acceptable
max_fraction_counts_excluded: 0.95
use_seed = False
num_iterations = 1000
# Set-up for pipeline notebooks; do not modify unless you are a power user!
notebook_basenames_list: Dual CRISPR 6-Scoring Preparation.ipynb,Dual CRISPR 7-Abundance Thresholds.ipynb,Dual CRISPR 8-Construct Scoring.ipynb

[test]
use_seed = True
num_iterations = 2
"""

        count_fps_or_dirs = "/my/counts/dir1,/my/counts_dir2"
        day_timepoints_str = "3,10,21"
        provided_output_dir = "/my/output_parent"
        is_test = False

        # count_fps_or_dirs string, day_timepoints_str, time_prefixes, and notebooks list should NOT be parsed to list
        # use_seed SHOULD be parsed to boolean
        # num_iterations SHOULD be parsed to int
        expected_output = {'machine_configuration': "laptop",
                           'count_fps_or_dirs': '/my/counts/dir1,/my/counts_dir2',
                           'day_timepoints_str': '3,10,21',
                           'processed_data_dir': "/my/output_parent",
                           'notebook_basenames_list': 'Dual CRISPR 6-Scoring Preparation.ipynb,Dual '
                            'CRISPR 7-Abundance Thresholds.ipynb,Dual CRISPR '
                             '8-Construct Scoring.ipynb',
                           'min_count_limit': 10,
                           'max_fraction_acceptable_spline_density_diff': 0.02,
                           'max_fraction_counts_excluded': 0.95,
                           'num_iterations': 1000,
                           'time_prefixes': 'T,D',
                           'use_seed': False}

        temp_config = tempfile.NamedTemporaryFile(mode="w")
        temp_config.write(config_str)
        temp_config.seek(0)

        real_output = ns_test._set_params(count_fps_or_dirs, day_timepoints_str, provided_output_dir, is_test,
                                          temp_config.name)
        self.assertEqual(expected_output, real_output)

    # end region
