# standard libraries
import os
import unittest
import warnings

import ccbbucsd.utilities.config_loader as ns_config

# test library
import run_mali_scoring as ns_test


class TestFunctions(unittest.TestCase):
    # no tests for _parse_cmd_line_args as it is so simple
    # no tests for main as it just chains together calls to other tested methods

    # region _set_params
    def test__set_params_test(self):
        config_str = """[DEFAULT]
machine_configuration = c4_2xlarge
keep_gzs: False
full_5p_r1: TATATATCTTGTGGAAAGGACGAAACACCG
full_5p_r2: CCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC
full_3p_r1: GTTTCAGAGCTATGCTGGAAACTGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAG
full_3p_r2: CAAACAAGGCTTTTCTCCAAGGGATATTTATAGTCTCAAAACACACAATTACTTTACAGTTAGGGTGAGTTTCCTTTTGTGCTGTTTTTTAAAATA
len_of_seq_to_match = 19
num_allowed_mismatches = 1
# min_count_limit in absolute counts, not log2
min_count_limit: 10
# max_fraction_acceptable_spline_density_diff is % of diff between max spline and min density
max_fraction_acceptable_spline_density_diff: 0.02
# any threshold throwing out > max_fraction_counts_excluded% of counts is not acceptable
max_fraction_counts_excluded: 0.95
use_seed = False
num_iterations = 1000
# Set-up for pipeline notebooks; do not modify unless you are a power user!
time_prefixes: T,D
count_notebooks: Dual CRISPR 1-Construct Scaffold Trimming.ipynb,Dual CRISPR 2-Constuct Filter.ipynb,Dual CRISPR 3-Construct Counting.ipynb,Dual CRISPR 4-Count Combination.ipynb,Dual CRISPR 5-Count Plots.ipynb
score_notebooks: Dual CRISPR 6-Scoring Preparation.ipynb,Dual CRISPR 7-Abundance Thresholds.ipynb,Dual CRISPR 8-Construct Scoring.ipynb

[c4_2xlarge]
main_dir: /home/ec2-user
data_dir: /data
num_processors: 7
# Set-up for directory structure; do not modify unless you are a power user!
code_dir: ${main_dir}/src/python
notebook_dir: ${main_dir}/notebooks
libraries_dir: ${main_dir}/library_definitions
raw_data_dir: ${data_dir}/raw
interim_data_dir: ${data_dir}/interim
processed_data_dir: ${data_dir}/processed

[laptop]
main_dir: /Users/Birmingham/Work/Repositories/ccbb_tickets_2017/
data_dir: /Users/Birmingham/Work/Data
num_processors: 3
# Set-up for directory structure; do not modify unless you are a power user!
code_dir: ${main_dir}/src/python
notebook_dir: ${main_dir}/notebooks
libraries_dir: ${main_dir}/library_definitions
raw_data_dir: ${data_dir}/raw
interim_data_dir: ${data_dir}/interim
processed_data_dir: ${data_dir}/processed

[test]
use_seed = True
num_iterations = 2
    """

        count_fps_or_dirs = "/my/counts/dir1,/my/counts_dir2"
        day_timepoints_str = "3,10,21"
        is_test = True

        # count_fps_or_dirs string, day_timepoints_str, time_prefixes, and notebooks list should NOT be parsed to list
        # use_seed SHOULD be parsed to boolean
        # num_iterations SHOULD be parsed to int
        expected_output = {'count_fps_or_dirs': '/my/counts/dir1,/my/counts_dir2',
                           'day_timepoints_str': '3,10,21',
                           'notebook_basenames_list': 'Dual CRISPR 6-Scoring Preparation.ipynb,Dual '
                                                      'CRISPR 7-Abundance Thresholds.ipynb,Dual CRISPR '
                                                      '8-Construct Scoring.ipynb',
                           'num_iterations': 2,
                           'time_prefixes': 'T,D',
                           'use_seed': True}

        temp_configfile_fp = ns_config.get_default_config_fp()
        try:
            with open(ns_config.get_default_config_fp(), "w") as f:
                f.write(config_str)

            with warnings.catch_warnings(record=True) as warnings_list:
                # Cause all warnings to always be triggered.
                warnings.simplefilter("always")

                real_output = ns_test._set_params(count_fps_or_dirs, day_timepoints_str, is_test)

                assert len(warnings_list) == 1
                assert "Scoring is running in TEST MODE; do not use results for data analysis!" in str(warnings_list[-1].message)


        finally:
            if os.path.exists(temp_configfile_fp):
                os.remove(temp_configfile_fp)

        self.assertEqual(expected_output, real_output)

    def test_set_params_real(self):
        config_str = """[DEFAULT]
machine_configuration = c4_2xlarge
keep_gzs: False
full_5p_r1: TATATATCTTGTGGAAAGGACGAAACACCG
full_5p_r2: CCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC
full_3p_r1: GTTTCAGAGCTATGCTGGAAACTGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAG
full_3p_r2: CAAACAAGGCTTTTCTCCAAGGGATATTTATAGTCTCAAAACACACAATTACTTTACAGTTAGGGTGAGTTTCCTTTTGTGCTGTTTTTTAAAATA
len_of_seq_to_match = 19
num_allowed_mismatches = 1
# min_count_limit in absolute counts, not log2
min_count_limit: 10
# max_fraction_acceptable_spline_density_diff is % of diff between max spline and min density
max_fraction_acceptable_spline_density_diff: 0.02
# any threshold throwing out > max_fraction_counts_excluded% of counts is not acceptable
max_fraction_counts_excluded: 0.95
use_seed = False
num_iterations = 1000
# Set-up for pipeline notebooks; do not modify unless you are a power user!
time_prefixes: T,D
count_notebooks: Dual CRISPR 1-Construct Scaffold Trimming.ipynb,Dual CRISPR 2-Constuct Filter.ipynb,Dual CRISPR 3-Construct Counting.ipynb,Dual CRISPR 4-Count Combination.ipynb,Dual CRISPR 5-Count Plots.ipynb
score_notebooks: Dual CRISPR 6-Scoring Preparation.ipynb,Dual CRISPR 7-Abundance Thresholds.ipynb,Dual CRISPR 8-Construct Scoring.ipynb

[c4_2xlarge]
main_dir: /home/ec2-user
data_dir: /data
num_processors: 7
# Set-up for directory structure; do not modify unless you are a power user!
code_dir: ${main_dir}/src/python
notebook_dir: ${main_dir}/notebooks
libraries_dir: ${main_dir}/library_definitions
raw_data_dir: ${data_dir}/raw
interim_data_dir: ${data_dir}/interim
processed_data_dir: ${data_dir}/processed

[laptop]
main_dir: /Users/Birmingham/Work/Repositories/ccbb_tickets_2017/
data_dir: /Users/Birmingham/Work/Data
num_processors: 3
# Set-up for directory structure; do not modify unless you are a power user!
code_dir: ${main_dir}/src/python
notebook_dir: ${main_dir}/notebooks
libraries_dir: ${main_dir}/library_definitions
raw_data_dir: ${data_dir}/raw
interim_data_dir: ${data_dir}/interim
processed_data_dir: ${data_dir}/processed

[test]
use_seed = True
num_iterations = 2
"""

        count_fps_or_dirs = "/my/counts/dir1,/my/counts_dir2"
        day_timepoints_str = "3,10,21"
        is_test = False

        # count_fps_or_dirs string, day_timepoints_str, time_prefixes, and notebooks list should NOT be parsed to list
        # use_seed SHOULD be parsed to boolean
        # num_iterations SHOULD be parsed to int
        expected_output = {'count_fps_or_dirs': '/my/counts/dir1,/my/counts_dir2',
                           'day_timepoints_str': '3,10,21',
                           'notebook_basenames_list': 'Dual CRISPR 6-Scoring Preparation.ipynb,Dual '
                            'CRISPR 7-Abundance Thresholds.ipynb,Dual CRISPR '
                             '8-Construct Scoring.ipynb',
                           'num_iterations': 1000,
                           'time_prefixes': 'T,D',
                           'use_seed': False}

        temp_configfile_fp = ns_config.get_default_config_fp()
        try:
            with open(ns_config.get_default_config_fp(), "w") as f:
                f.write(config_str)

            real_output = ns_test._set_params(count_fps_or_dirs, day_timepoints_str, is_test)
        finally:
            if os.path.exists(temp_configfile_fp):
                os.remove(temp_configfile_fp)

        self.assertEqual(expected_output, real_output)

    # end region
