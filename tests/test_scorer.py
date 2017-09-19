# standard libraries
import unittest

# third-party libraries
import pandas
import pandas.util.testing

# project-specific libraries
import dual_crispr.scorer as ns_test
import tests.test_fitted_per_replicate_data as ns_fitted
import tests.test_thresholded_per_replicate_data as ns_thresh
import tests.test_per_replicate_data as ns_help

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
    # def test_get_constructs_to_ignore_across_replicates_mask(self):
    #     self.fail("test not implemented")

    def test_generate_naive_fitness_per_construct_series(self):
        input_thresh_rep_list = [ns_thresh.TestThresholdedPerReplicateData.help_make_thresh_rep(True),
                                 ns_thresh.TestThresholdedPerReplicateData.help_make_thresh_rep(False)]
        expected_output_series = ns_help.TestScorer.help_make_naive_fcs_series()

        real_output_series = ns_test._generate_naive_fitness_per_construct_series(input_thresh_rep_list)
        self.assertEqual(len(expected_output_series), len(real_output_series))
        error_msgs = []

        for i in range(len(expected_output_series)):
            if abs(expected_output_series[i] - real_output_series[i]) > 0.00051:
                error_msgs.append(
                    "For probe pair id '{0}', expected value {1} != real value {2} within 0.0005 delta".format(
                        expected_output_series.index[i], expected_output_series[i], real_output_series[i]))

        print("\n".join(error_msgs))
        self.assertEqual(len(error_msgs), 0)

    def test__make_multiindex_df_across_replicates(self):
        input_thresh_rep_list = [ns_thresh.TestThresholdedPerReplicateData.help_make_thresh_rep(True),
                                 ns_thresh.TestThresholdedPerReplicateData.help_make_thresh_rep(False)]
        expected_output_df = ns_help.TestScorer._help_make_multiindex_df()
        real_output_df = ns_test._make_multiindex_df_across_replicates(input_thresh_rep_list, "descriptive_stats_df")

        # belt and suspenders approach: check values within a precision with code below, then check indexes, etc,
        # with assert_frame_equal below (with unrealistically strong rounding to make values match)
        ns_help.help_test_df_equality(self, expected_output_df, real_output_df, 0.0000005)

    def test__calc_sum_of_stat_per_construct_across_replicates_series(self):
        input_df = ns_help.TestScorer._help_make_multiindex_df()
        expected_output_series = ns_help.access_test_file_as_series("sum_covariance_of_log2_fractions_w_timepoints.txt",
                                                                    set_index=True)
        real_output_series = ns_test._calc_sum_of_stat_per_construct_across_replicates_series(
            input_df, "covariance_of_log2_fractions_w_timepoints")

        # belt and suspenders approach: check values within a precision with code below, then check indexes, etc,
        # with assert_frame_equal below (with unrealistically strong rounding to make values match)
        ns_help.help_test_series_equality(self, expected_output_series, real_output_series, 0.0000000005)

    def test__generate_stddev_of_naive_fitness_per_construct_series(self):
        fitted_per_rep_list = [ns_fitted.TestFittedPerReplicateData.help_make_fitted_rep(get_rep_1=True),
                               ns_fitted.TestFittedPerReplicateData.help_make_fitted_rep(get_rep_1=False)]
        ns_test._generate_stddev_of_naive_fitness_per_construct_series(fitted_per_rep_list)
        self.fail("test not implemented")
