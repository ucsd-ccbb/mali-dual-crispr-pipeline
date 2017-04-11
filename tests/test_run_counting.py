# standard libraries
import tempfile
import unittest

from dual_crispr import run_counting as ns_test


class TestFunctions(unittest.TestCase):
    # no tests for _parse_cmd_line_args as it is so simple
    # no tests for main as it just chains together calls to other tested methods

    def test__set_params(self):
        config_str = """[DEFAULT]
machine_configuration = c4_2xlarge

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
        input_fastq_dir_name = "/data/raw/test_fastq_dir"
        expected_output_dir_path = '/data/output'

        expected_output = {'machine_configuration': 'c4_2xlarge',
                           'raw_data_dir': '/data/raw/test_fastq_dir',
                           'interim_data_dir': '/data/output/temporary_files/{run_prefix}',
                           'processed_data_dir': '/data/output',
                           'notebook_basenames_list': 'Dual CRISPR 1-Construct Scaffold Trimming.ipynb,Dual CRISPR '
                                                      '2-Constuct Filter.ipynb,Dual CRISPR 3-Construct '
                                                      'Counting.ipynb,Dual CRISPR 4-Count Combination.ipynb,'
                                                      'Dual CRISPR 5-Count Plots.ipynb',
                           'full_5p_r1': 'TATATATCTTGTGGAAAGGACGAAACACCG',
                           'full_5p_r2': 'CCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC',
                           'full_3p_r1': 'GTTTCAGAGCTATGCTGGAAACTGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAG',
                           'full_3p_r2': 'CAAACAAGGCTTTTCTCCAAGGGATATTTATAGTCTCAAAACACACAATTACTTTACAGTTAGGGTGAGTTTCCTTTTGTGCTGTTTTTTAAAATA',
                           'len_of_seq_to_match': 19,
                           'num_allowed_mismatches': 1,
                           'g_fastqs_dir':'/data/raw/test_fastq_dir',
                           'g_trimmed_fastqs_dir':'/data/output/temporary_files/{run_prefix}',
                           'g_filtered_fastqs_dir': '/data/output/temporary_files/{run_prefix}'
                           }

        temp_config = tempfile.NamedTemporaryFile(mode="w")
        temp_config.write(config_str)
        temp_config.seek(0)

        real_output = ns_test._set_params(input_fastq_dir_name, expected_output_dir_path, temp_config.name)
        self.maxDiff = None
        self.assertEqual(expected_output, real_output)