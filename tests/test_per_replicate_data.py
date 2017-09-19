# standard libraries
import io
import os
import unittest

# third-party libraries
import pandas
import pandas.util.testing

# project-specific libraries
import dual_crispr.per_replicate_data as ns_test

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


def access_test_file_as_df(file_name, sub_folder=None, col_index_to_use_as_index=0):
    """

    Args:
        file_name:
        sub_folder:
        col_index_to_use_as_index (:obj:`int`, optional):

    Returns:
        pandas.DataFrame
    """
    test_data_fp = get_test_file_path(file_name, sub_folder)
    test_data_df = pandas.read_csv(test_data_fp, sep="\t", index_col=col_index_to_use_as_index)
    return test_data_df


def access_test_file_as_series(file_name, sub_folder=None, clear_index_name=False, set_index=False):
    test_data_df = access_test_file_as_df(file_name, sub_folder=sub_folder, col_index_to_use_as_index=0)
    return _make_series_from_df(test_data_df, clear_index_name=clear_index_name, set_index=set_index)


def get_test_file_path(file_name, sub_folder=None):
    if sub_folder is None: sub_folder = get_set_8_dir_name()
    dir_of_test_script = os.path.dirname(os.path.realpath(__file__))
    test_data_fp = os.path.join(dir_of_test_script, "../dual_crispr/distributed_files/test_data",
                                sub_folder, file_name)
    return test_data_fp


def help_make_df(contents_str, col_index_to_use_as_index=None):
    """

    Args:
        contents_str:
        col_index_to_use_as_index:

    Returns:
        pandas.DataFrame
    """
    df_file_obj = io.StringIO(contents_str)
    return pandas.read_csv(df_file_obj, sep="\t", index_col=col_index_to_use_as_index)


def help_make_series(contents_str, col_index_to_use_as_index=None, clear_index_name=False, set_index=False):
    test_data_df = help_make_df(contents_str, col_index_to_use_as_index=col_index_to_use_as_index)
    return _make_series_from_df(test_data_df, clear_index_name=clear_index_name, set_index=set_index)


def _make_series_from_df(test_data_df, clear_index_name=False, set_index=False):
    test_data_series = test_data_df.ix[:, 0]
    if set_index:
        test_data_series.index = test_data_df.index

    test_data_series.name = None
    if clear_index_name:
        test_data_series.index.names = [None]
    return test_data_series


def get_set_8_dir_name():
    return "test_set_8"


def help_test_df_equality(test_case, expected_output_df, real_output_df, max_allowable_diff_per_field, round_to=-10):
    # belt and suspenders approach: check values within a precision with code below, then check indexes, etc,
    # with assert_frame_equal below (with unrealistically strong rounding to make values always match)
    differences_df = pandas.DataFrame.abs(expected_output_df - real_output_df) > max_allowable_diff_per_field
    # add all trues in dataframe on 0 axis to get series, then add all numbers in series to get total
    # should come out to be zero if no differences beyond the precision
    num_different = differences_df.sum(axis=0).sum()
    test_case.assertEqual(0, num_different)
    pandas.util.testing.assert_frame_equal(expected_output_df.round(decimals=round_to),
                                           real_output_df.round(decimals=round_to))


def help_test_series_equality(test_case, expected_output_series, real_output_series, max_allowable_diff_per_field,
                              round_to=-10):
    # belt and suspenders approach: check values within a precision with code below, then check indexes, etc,
    # with assert_frame_equal below (with unrealistically strong rounding to make values always match)
    differences_series = pandas.Series.abs(expected_output_series - real_output_series) > max_allowable_diff_per_field
    # add all trues in seies to get total
    # should come out to be zero if no differences beyond the precision
    num_different = differences_series.sum()
    test_case.assertEqual(0, num_different)
    pandas.util.testing.assert_series_equal(expected_output_series.round(decimals=round_to),
                                            real_output_series.round(decimals=round_to))


class TestPerReplicateData(unittest.TestCase):
    @staticmethod
    def help_make_timepoints_series(get_rep_1=True):
        if get_rep_1:
            result =  pandas.Series(data=[21, 28], index=["A549CV4_T21_1", "A549CV4_T28_1"])
        else:
            raise NotImplementedError("not implemented for replicate 2")

        return result

    @staticmethod
    def help_make_per_rep_data(get_rep_1=True, input_log2_fractions_df=None):
        """

        Args:
            get_rep_1:
            input_log2_fractions_df:

        Returns:
            PerReplicateData
        """
        if input_log2_fractions_df is None:
            input_log2_fractions_df = \
                TestPerReplicateData.help_make_log2_fractions_per_construct_per_timepoint_df(get_rep_1)

        if get_rep_1:
            per_rep_data = ns_test.PerReplicateData(1, input_log2_fractions_df, pandas.Series(
                data=[-18.45836, -19.32780], index=pandas.Index(
                    data=["A549CV4_T21_1", "A549CV4_T28_1"], name="sampleName")))
        else:
            per_rep_data = ns_test.PerReplicateData(2, input_log2_fractions_df, pandas.Series(
                data=[-18.66791, -19.97675], index=pandas.Index(
                    data=["A549CV4_T21_2", "A549CV4_T28_2"], name="sampleName")))

        return per_rep_data

    @staticmethod
    def help_make_log2_fractions_per_construct_per_timepoint_df(get_rep_1=True):
        file_name = "5_postbelowthreshold_x1_new.txt" if get_rep_1 else "5_postbelowthreshold_x2_new.txt"
        return access_test_file_as_df(file_name)

    def test_extract_values_set_from_data_headers_replicates(self):
        input_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	153	321
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""
        input_df = help_make_df(input_str)
        real_output = ns_test.PerReplicateData.extract_values_set_from_data_headers(input_df, get_timepts=False)

        expected_output = [1, 2]
        self.assertListEqual(expected_output, real_output)

    def test_extract_values_set_from_data_headers_timepts(self):
        input_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	153	321
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""
        input_df = help_make_df(input_str)
        real_output = ns_test.PerReplicateData.extract_values_set_from_data_headers(input_df, get_timepts=True)

        expected_output = [21, 28]
        self.assertListEqual(expected_output, real_output)

    def test__timepoints_series(self):
        per_rep_data = self.help_make_per_rep_data()
        expected_output = pandas.Series(data=[21, 28], index=["A549CV4_T21_1", "A549CV4_T28_1"])
        pandas.util.testing.assert_series_equal(expected_output, per_rep_data.timepoints_series)


# Ok, yes, this class should be in test_scorer.py.  But it can't be because that causes a circular reference
# (since test_scorer has to reference test_fitted_per_replicate_data, but test_fitted_per_replicate_data has
# to ALSO reference test_scorer if test_scorer holds help_make_naive_fcs_series).
class TestScorer:
    @staticmethod
    def help_make_naive_fcs_series():
        """

        Returns:
            pandas.Series
        """
        return access_test_file_as_series("6_postposteriorprobplot_fc_new.txt", clear_index_name=True)

    @staticmethod
    def _help_make_multiindex_df():
        expected_output_fp = get_test_file_path("multiindex_descriptive_stats_corrected.txt")
        expected_output_df = pandas.read_csv(expected_output_fp, sep="\t", skiprows=1, index_col=0)
        # the multiindex in the test file contains the correct info, but pandas reads the integers levels in as
        # strings so the indices not actually match between expected and real, so manually rewrite the index
        expected_output_df.columns = pandas.MultiIndex(levels=[[1, 2], ['mean_usable_log2_fraction_for_construct',
                                                                        'mean_of_usable_timepoints',
                                                                        'variance_of_usable_timepoints',
                                                                        'covariance_of_log2_fractions_w_timepoints']],
                                                       labels=[[0, 0, 0, 0, 1, 1, 1, 1], [0, 1, 2, 3, 0, 1, 2, 3]])
        return expected_output_df
