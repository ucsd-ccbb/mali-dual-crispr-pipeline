# standard libraries
import unittest

# third-party libraries
import pandas
import pandas.util.testing

# project-specific libraries
import dual_crispr.thresholded_per_replicate_data as ns_test
import tests.test_per_replicate_data as ns_help

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestThresholdedPerReplicateData(unittest.TestCase):
    @staticmethod
    def help_make_thresh_rep(get_rep_1=True):
        per_rep_data = ns_help.TestPerReplicateData.help_make_per_rep_data(get_rep_1)
        fail_mask_series, usable_constructs_mask_df = \
            TestThresholdedPerReplicateData._help_make_fail_mask_and_usable_constructs(get_rep_1)
        descriptive_stats_df = TestThresholdedPerReplicateData._help_make_descriptive_stats_df(get_rep_1)

        result = ns_test.ThresholdedPerReplicateData(per_rep_data.replicate_num,
                                                     per_rep_data._log2_fractions_by_constructs_by_timepoints_df,
                                                     per_rep_data._log2_fractions_thresholds_by_timepoints_series, 2,
                                                     _run_init_methods=False)
        result._fails_min_num_timepts_bool_by_construct_series = fail_mask_series
        result._usable_constructs_by_timepoint_mask = usable_constructs_mask_df
        result._descriptive_stats_by_construct_df = descriptive_stats_df
        return result

    @staticmethod
    def _help_make_fail_mask_series(get_rep_1=True):
        file_name = "5_postbelowthreshold_bad1_new.txt" if get_rep_1 else "5_postbelowthreshold_bad2_new.txt"
        fail_mask_series = ns_help.access_test_file_as_series(file_name)
        return fail_mask_series

    @staticmethod
    def _help_make_usable_constructs_by_timepoint_mask_df(get_rep_1=True):
        file_name = "5_postbelowthreshold_good1_new.txt" if get_rep_1 else "5_postbelowthreshold_good2_new.txt"
        return ns_help.access_test_file_as_df(file_name)

    @staticmethod
    def _help_make_fail_mask_and_usable_constructs(get_rep_1=True):
        fail_mask_series = TestThresholdedPerReplicateData._help_make_fail_mask_series(get_rep_1)
        usable_constructs_by_timepoint_mask_df = \
            TestThresholdedPerReplicateData._help_make_usable_constructs_by_timepoint_mask_df(get_rep_1)
        return fail_mask_series, usable_constructs_by_timepoint_mask_df

    @staticmethod
    def _help_make_descriptive_stats_df(get_rep_1=True):
        stats_file_name = "rep1_descriptive_stats.txt" if get_rep_1 else "rep2_descriptive_stats.txt"
        descriptive_stats_df = ns_help.access_test_file_as_df(stats_file_name)
        return descriptive_stats_df


class TestFunctions(unittest.TestCase):
    # No test for generate_thresholded_per_rep_data_from_per_rep_data because it just passes through values to init

    def test__get_passes_threshold_bool_by_construct_by_timept_df(self):
        input_log2_fractions_by_constructs_by_timepoints_str = """probe_pair_id	A549CV4_T21_1	A549CV4_T28_1
0NonTargetingControlGuideForHuman0412_BRCA1_chr17_41276018	-22.98336	-23.252803
0NonTargetingControlGuideForHuman0352_SETD2_chr3_47142972	-12.77757	-14.136459
0NonTargetingControlGuideForHuman0362_BRCA1_chr17_41256141	-15.39091	-15.995415
0NonTargetingControlGuideForHuman0352_BRCA2_chr13_32912006	-19.05263	-16.608947
0NonTargetingControlGuideForHuman0362_RB1_chr13_48934148	-15.49151	-16.082878
0NonTargetingControlGuideForHuman0412_FLT3_chr13_28636131	-15.49955	-16.175988
0NonTargetingControlGuideForHuman0352_BRCA1_chr17_41256141	-13.14889	-21.339914
0NonTargetingControlGuideForHuman0362_CDKN2A_chr9_21968726	-15.34674	-15.601752
0NonTargetingControlGuideForHuman0412_HDAC1_chr1_32757816	-15.96099	-16.310289
0NonTargetingControlGuideForHuman0362_SMO_chr7_128845117	-16.75454	-16.949023
0NonTargetingControlGuideForHuman0352_AKT1_chr14_105239686	-13.65918	-13.37322
"""
        expected_output_str = """probe_pair_id	A549CV4_T21_1	A549CV4_T28_1
0NonTargetingControlGuideForHuman0412_BRCA1_chr17_41276018	False	False
0NonTargetingControlGuideForHuman0352_SETD2_chr3_47142972	True	True
0NonTargetingControlGuideForHuman0362_BRCA1_chr17_41256141	True	True
0NonTargetingControlGuideForHuman0352_BRCA2_chr13_32912006	False	True
0NonTargetingControlGuideForHuman0362_RB1_chr13_48934148	True	True
0NonTargetingControlGuideForHuman0412_FLT3_chr13_28636131	True	True
0NonTargetingControlGuideForHuman0352_BRCA1_chr17_41256141	True	False
0NonTargetingControlGuideForHuman0362_CDKN2A_chr9_21968726	True	True
0NonTargetingControlGuideForHuman0412_HDAC1_chr1_32757816	True	True
0NonTargetingControlGuideForHuman0362_SMO_chr7_128845117	True	True
0NonTargetingControlGuideForHuman0352_AKT1_chr14_105239686	True	True
"""
        input_log2_fractions_by_constructs_by_timepoints_df = ns_help.help_make_df(
            input_log2_fractions_by_constructs_by_timepoints_str, col_index_to_use_as_index=0)
        expected_output_df = ns_help.help_make_df(expected_output_str, col_index_to_use_as_index=0)
        per_rep_data = ns_help.TestPerReplicateData.help_make_per_rep_data()

        real_output_df = ns_test._get_passes_threshold_bool_by_construct_by_timept_df(
            input_log2_fractions_by_constructs_by_timepoints_df,
            per_rep_data._log2_fractions_thresholds_by_timepoints_series)

        pandas.util.testing.assert_frame_equal(expected_output_df, real_output_df)

    def test__get_constructs_to_ignore_mask_series(self):
        expected_output_series = TestThresholdedPerReplicateData._help_make_fail_mask_series()
        per_rep_data = ns_help.TestPerReplicateData.help_make_per_rep_data()

        real_output_series = ns_test._get_constructs_to_ignore_mask_series(
            per_rep_data._log2_fractions_by_constructs_by_timepoints_df,
            per_rep_data._log2_fractions_thresholds_by_timepoints_series, 2)
        self.assertTrue(expected_output_series.equals(real_output_series))

    def test__get_usable_constructs_by_timepoints_mask_df(self):
        input_fail_mask_series, expected_output_df = \
            TestThresholdedPerReplicateData._help_make_fail_mask_and_usable_constructs()
        per_rep_data = ns_help.TestPerReplicateData.help_make_per_rep_data()

        real_output_df = ns_test._get_usable_constructs_by_timepoints_mask_df(
            per_rep_data._log2_fractions_by_constructs_by_timepoints_df,
            per_rep_data._log2_fractions_thresholds_by_timepoints_series,
            input_fail_mask_series)
        pandas.util.testing.assert_frame_equal(expected_output_df, real_output_df)

    def test__calculate_covariance_by_construct_series(self):
        input_usable_log2_fractions_by_constructs_str = """probe_pair_id	A549CV4_T21_1	A549CV4_T28_1
0NonTargetingControlGuideForHuman0412_BRCA1_chr17_41276018	NaN	NaN
0NonTargetingControlGuideForHuman0352_SETD2_chr3_47142972	-12.77757	-14.136459
0NonTargetingControlGuideForHuman0362_BRCA1_chr17_41256141	-15.39091	-15.995415
0NonTargetingControlGuideForHuman0352_BRCA2_chr13_32912006	NaN	NaN
0NonTargetingControlGuideForHuman0362_RB1_chr13_48934148	-15.49151	-16.082878
0NonTargetingControlGuideForHuman0412_FLT3_chr13_28636131	-15.49955	-16.175988
0NonTargetingControlGuideForHuman0352_BRCA1_chr17_41256141	NaN	NaN
0NonTargetingControlGuideForHuman0362_CDKN2A_chr9_21968726	-15.34674	-15.601752
0NonTargetingControlGuideForHuman0412_HDAC1_chr1_32757816	-15.96099	-16.310289
0NonTargetingControlGuideForHuman0362_SMO_chr7_128845117	-16.75454	-16.949023
0NonTargetingControlGuideForHuman0352_AKT1_chr14_105239686	-13.65918	-13.37322
"""
        input_usable_timepoints_by_construct_str = """probe_pair_id	A549CV4_T21_1	A549CV4_T28_1
0NonTargetingControlGuideForHuman0412_BRCA1_chr17_41276018	NaN	NaN
0NonTargetingControlGuideForHuman0352_SETD2_chr3_47142972	21	28
0NonTargetingControlGuideForHuman0362_BRCA1_chr17_41256141	21	28
0NonTargetingControlGuideForHuman0352_BRCA2_chr13_32912006	NaN	NaN
0NonTargetingControlGuideForHuman0362_RB1_chr13_48934148	21	28
0NonTargetingControlGuideForHuman0412_FLT3_chr13_28636131	21	28
0NonTargetingControlGuideForHuman0352_BRCA1_chr17_41256141	NaN	NaN
0NonTargetingControlGuideForHuman0362_CDKN2A_chr9_21968726	21	28
0NonTargetingControlGuideForHuman0412_HDAC1_chr1_32757816	21	28
0NonTargetingControlGuideForHuman0362_SMO_chr7_128845117	21	28
0NonTargetingControlGuideForHuman0352_AKT1_chr14_105239686	21	28
"""
        expected_output_str = """probe_pair_id	
0NonTargetingControlGuideForHuman0412_BRCA1_chr17_41276018	0
0NonTargetingControlGuideForHuman0352_SETD2_chr3_47142972	-2.378057095
0NonTargetingControlGuideForHuman0362_BRCA1_chr17_41256141	-1.057891932
0NonTargetingControlGuideForHuman0352_BRCA2_chr13_32912006	0
0NonTargetingControlGuideForHuman0362_RB1_chr13_48934148	-1.034895007
0NonTargetingControlGuideForHuman0412_FLT3_chr13_28636131	-1.183771156
0NonTargetingControlGuideForHuman0352_BRCA1_chr17_41256141	0
0NonTargetingControlGuideForHuman0362_CDKN2A_chr9_21968726	-0.446273467
0NonTargetingControlGuideForHuman0412_HDAC1_chr1_32757816	-0.611264129
0NonTargetingControlGuideForHuman0362_SMO_chr7_128845117	-0.34033724
0NonTargetingControlGuideForHuman0352_AKT1_chr14_105239686	0.500433889
"""
        input_usable_log2_fractions_by_constructs_df = ns_help.help_make_df(
            input_usable_log2_fractions_by_constructs_str, col_index_to_use_as_index=0)
        input_usable_timepoints_by_construct_df = ns_help.help_make_df(
            input_usable_timepoints_by_construct_str, col_index_to_use_as_index=0)
        expected_output_series = ns_help.help_make_series(expected_output_str, col_index_to_use_as_index=0)

        real_output_series = ns_test._calculate_covariance_by_construct_series(
            input_usable_log2_fractions_by_constructs_df, input_usable_timepoints_by_construct_df)

        pandas.util.testing.assert_series_equal(expected_output_series, real_output_series, check_less_precise=True)

    def test__generate_descriptive_statistics_by_construct(self):
        input_fail_mask_series, input_usable_constructs_by_timepoint_mask_df = \
            TestThresholdedPerReplicateData._help_make_fail_mask_and_usable_constructs()
        per_rep_data = ns_help.TestPerReplicateData.help_make_per_rep_data()
        expected_output_df = TestThresholdedPerReplicateData._help_make_descriptive_stats_df()

        real_output_df = ns_test._generate_descriptive_statistics_by_construct(
            per_rep_data._log2_fractions_by_constructs_by_timepoints_df,
            input_usable_constructs_by_timepoint_mask_df)

        pandas.util.testing.assert_frame_equal(expected_output_df, real_output_df)
