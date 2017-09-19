# standard libraries
import unittest

# third-party libraries
import pandas.util.testing

# project-specific libraries
import dual_crispr.fitted_per_replicate_data as ns_test
import tests.test_per_replicate_data as ns_help
import tests.test_thresholded_per_replicate_data as ns_thresh_help

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFittedPerReplicateData(unittest.TestCase):
    @staticmethod
    def help_make_fitted_rep(get_rep_1=True):
        thresh_per_rep_data = ns_thresh_help.TestThresholdedPerReplicateData.help_make_thresh_rep(get_rep_1=get_rep_1)
        naive_fitness_per_construct_series = ns_help.TestScorer.help_make_naive_fcs_series()

        result = ns_test.FittedPerReplicateData(thresh_per_rep_data.replicate_num,
                                                thresh_per_rep_data._log2_fractions_by_constructs_by_timepoints_df,
                                                thresh_per_rep_data._log2_fractions_thresholds_by_timepoints_series,
                                                thresh_per_rep_data._min_num_timepoints_above_abundance_threshold,
                                                naive_fitness_per_construct_series,
                                                _run_init_methods=False)

        result._normed_init_log2_fraction_per_construct_series = \
            TestFittedPerReplicateData.help_make_normed_init_log2_fraction_per_construct_series(get_rep_1)
        result._lambda_per_timepoint_series = \
            TestFittedPerReplicateData.help_make_lambda_by_timepoints_series(get_rep_1)
        result._expected_log2_fraction_per_construct_per_timepoint_df = \
            TestFittedPerReplicateData.help_make_expected_log2_fraction_per_construct_per_timepoint_df(get_rep_1)
        return result

    @staticmethod
    def help_make_normed_init_log2_fraction_per_construct_series(get_rep_1=True):
        file_name = "6_postposteriorprobplot_a1_new.txt" if get_rep_1 else "6_postposteriorprobplot_a2_new.txt"
        return ns_help.access_test_file_as_series(file_name)

    @staticmethod
    def help_make_lambda_by_timepoints_series(get_rep_1=True):
        if get_rep_1:
            file_name = "lambda1.txt"  # "lambda_of_time_mult_by_6_postposteriorprobplot_fc_new_plus_6_postposteriorprobplot_a1_new.txt"
        else:
            file_name = "lambda2.txt"  # "lambda_of_time_mult_by_6_postposteriorprobplot_fc_new_plus_6_postposteriorprobplot_a2_new.txt"
        return ns_help.access_test_file_as_series(file_name, clear_index_name=True)

    @staticmethod
    def help_make_expected_log2_fraction_per_construct_per_timepoint_df(get_rep_1=True):
        file_name = "xfit1.txt" if get_rep_1 else "xfit2.txt"
        return ns_help.access_test_file_as_df(file_name)


class TestFunctions(unittest.TestCase):
    def test__calculate_normed_init_log2_fraction_per_construct_series(self):
        # only tested for replicate 1
        input_descriptive_stats_df = ns_help.access_test_file_as_df("rep1_descriptive_stats.txt")
        input_naive_fcs_series = ns_help.TestScorer.help_make_naive_fcs_series()
        input_log2_fractions_by_constructs_by_tmpts_df = \
            ns_help.TestPerReplicateData.help_make_log2_fractions_per_construct_per_timepoint_df()
        input_failed_constructs_series = \
            ns_thresh_help.TestThresholdedPerReplicateData._help_make_fail_mask_series()

        expected_output_series = TestFittedPerReplicateData.help_make_normed_init_log2_fraction_per_construct_series()

        real_output_series = ns_test._calculate_normed_init_log2_fraction_per_construct_series(
            input_descriptive_stats_df, input_naive_fcs_series, input_log2_fractions_by_constructs_by_tmpts_df,
            input_failed_constructs_series)

        ns_help.help_test_series_equality(self, expected_output_series, real_output_series, 0.0125)

    def test__calculate_naive_fitness_mult_by_time_per_construct_per_timept_df(self):
        # only tested for replicate 1
        input_naive_fcs_series = ns_help.TestScorer.help_make_naive_fcs_series()
        input_timepoints_series = ns_help.TestPerReplicateData.help_make_timepoints_series()
        expected_output_df = ns_help.access_test_file_as_df("time_mult_by_6_postposteriorprobplot_fc_new.txt")

        real_output_df = ns_test._calculate_naive_fitness_mult_by_time_per_construct_per_timept_df(
            input_naive_fcs_series, input_timepoints_series)
        ns_help.help_test_df_equality(self, expected_output_df, real_output_df, 0.000001)

    def test__calculate_normed_init_log2_fraction_plus_naive_fitness_mult_by_time_per_construct_per_timept_df(self):
        # only tested for replicate 1
        input_normed_init_log2_fraction_per_construct_series = \
            TestFittedPerReplicateData.help_make_normed_init_log2_fraction_per_construct_series()
        input_naive_fcs_series = ns_help.TestScorer.help_make_naive_fcs_series()
        input_timepoints_series = ns_help.TestPerReplicateData.help_make_timepoints_series()
        expected_output_df = ns_help.access_test_file_as_df(
            "time_mult_by_6_postposteriorprobplot_fc_new_plus_6_postposteriorprobplot_a1_new.txt")

        real_output_df = \
            ns_test._calculate_normed_init_log2_fraction_plus_naive_fitness_mult_by_time_per_construct_per_timept_df(
                input_normed_init_log2_fraction_per_construct_series, input_naive_fcs_series, input_timepoints_series)
        pandas.util.testing.assert_frame_equal(expected_output_df, real_output_df)

    def test__calculate_lambda_per_timepoint_series(self):
        input_normed_init_log2_fraction_per_construct_series = \
            TestFittedPerReplicateData.help_make_normed_init_log2_fraction_per_construct_series()
        input_naive_fcs_series = ns_help.TestScorer.help_make_naive_fcs_series()
        input_timepoints_series = ns_help.TestPerReplicateData.help_make_timepoints_series()

        expected_output_series = TestFittedPerReplicateData.help_make_lambda_by_timepoints_series()
        real_output_series = ns_test._calculate_lambda_per_timepoint_series(
            input_normed_init_log2_fraction_per_construct_series, input_naive_fcs_series, input_timepoints_series)

        ns_help.help_test_series_equality(self, expected_output_series, real_output_series, 0.0002)
        # pandas.util.testing.assert_series_equal(expected_output_series, real_output_series)

    def test__calculate_expected_log2_fraction_per_construct_per_timepoint_df(self):
        input_timepoints_series = ns_help.TestPerReplicateData.help_make_timepoints_series()
        input_normed_init_log2_fraction_per_construct_series = \
            TestFittedPerReplicateData.help_make_normed_init_log2_fraction_per_construct_series()
        input_naive_fcs_series = ns_help.TestScorer.help_make_naive_fcs_series()
        input_lambda_by_timepoints_series = TestFittedPerReplicateData.help_make_lambda_by_timepoints_series()

        expected_output_df = \
            TestFittedPerReplicateData.help_make_expected_log2_fraction_per_construct_per_timepoint_df()

        real_output_df = ns_test._calculate_expected_log2_fraction_per_construct_per_timepoint_df(
            input_timepoints_series, input_normed_init_log2_fraction_per_construct_series, input_naive_fcs_series,
            input_lambda_by_timepoints_series)

        ns_help.help_test_df_equality(self, expected_output_df, real_output_df, 0.0145)

    def test_generate_fitted_per_rep_data_from_thresholded_per_rep_data(self):
        input_thresholded_rep = ns_thresh_help.TestThresholdedPerReplicateData.help_make_thresh_rep(get_rep_1=True)
        input_naive_fcs_series = ns_help.TestScorer.help_make_naive_fcs_series()

        expected_output_init_log2_fraction_per_construct_series = \
            TestFittedPerReplicateData.help_make_normed_init_log2_fraction_per_construct_series()
        expected_output_lambda_by_timepoints_series = \
            TestFittedPerReplicateData.help_make_lambda_by_timepoints_series()
        expected_output_expected_log2_fraction_per_construct_per_timepoint_df = \
            TestFittedPerReplicateData.help_make_expected_log2_fraction_per_construct_per_timepoint_df()

        output_obj = ns_test.generate_fitted_per_rep_data_from_thresholded_per_rep_data(input_thresholded_rep,
                                                                                             input_naive_fcs_series)

        ns_help.help_test_series_equality(self, expected_output_init_log2_fraction_per_construct_series,
                                          output_obj._normed_init_log2_fraction_per_construct_series, 0.0125)
        ns_help.help_test_series_equality(self, expected_output_lambda_by_timepoints_series,
                                                output_obj._lambda_per_timepoint_series, 0.00035)
        ns_help.help_test_df_equality(self, expected_output_expected_log2_fraction_per_construct_per_timepoint_df,
                                      output_obj._expected_log2_fraction_per_construct_per_timepoint_df, 0.0145)
