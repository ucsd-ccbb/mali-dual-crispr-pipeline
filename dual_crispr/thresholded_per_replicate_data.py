# third-party libraries
import numpy
import pandas

# ccbb libraries
import dual_crispr.scoring_prep as ns_prep
from dual_crispr.per_replicate_data import PerReplicateData

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


# NB: _run_init_methods param only exists to facilitate unit-testing and should NOT be changed for any other usage
class ThresholdedPerReplicateData(PerReplicateData):
    def __init__(self, replicate_num, log2_fractions_by_constructs_by_timepoints_df,
                 log2_fractions_thresholds_by_timepoints_series, min_num_timepoints_above_abundance_threshold,
                 _run_init_methods=True):
        super().__init__(replicate_num, log2_fractions_by_constructs_by_timepoints_df,
                         log2_fractions_thresholds_by_timepoints_series)
        self._min_num_timepoints_above_abundance_threshold = min_num_timepoints_above_abundance_threshold

        if _run_init_methods:
            self._fails_min_num_timepts_bool_by_construct_series = _get_constructs_to_ignore_mask_series(
                self._log2_fractions_by_constructs_by_timepoints_df, self._log2_fractions_thresholds_by_timepoints_series,
                self._min_num_timepoints_above_abundance_threshold)
            self._usable_constructs_by_timepoint_mask = _get_usable_constructs_by_timepoints_mask_df(
                self._log2_fractions_by_constructs_by_timepoints_df, self._log2_fractions_thresholds_by_timepoints_series,
                self._fails_min_num_timepts_bool_by_construct_series)
            self._descriptive_stats_by_construct_df = _generate_descriptive_statistics_by_construct(
                self._log2_fractions_by_constructs_by_timepoints_df, self._usable_constructs_by_timepoint_mask)

    @property
    def descriptive_stats_df(self):
        return self._descriptive_stats_by_construct_df.copy()


# Not a class method of ThresholdedPerReplicateData because anything that inherits from ThresholdedPerReplicateData
# shouldn't have this method
def generate_thresholded_per_rep_data_from_per_rep_data(
        per_replicate_data, min_num_timepoints_above_abundance_threshold):
    """

    Args:
        min_num_timepoints_above_abundance_threshold (integer):
        per_replicate_data (PerReplicateData):
    """
    result = ThresholdedPerReplicateData(per_replicate_data.replicate_num,
                                         per_replicate_data._log2_fractions_by_constructs_by_timepoints_df,
                                         per_replicate_data._log2_fractions_thresholds_by_timepoints_series,
                                         min_num_timepoints_above_abundance_threshold)
    return result


# Note: the below functions are NOT methods of the ThresholdedRpeplicateData object--even though that is the only
# entity that uses their functionality--because they should not be available to a ThresholdedReplicateData object
# after it is created.  That is, once you've run _generate_descriptive_statistics_by_construct, you should not be
# running it again on the ThresholdedPerReplicateData object later, and more importantly, a FittedPerReplicateData
# object derived from a ThresholdedPerReplicateData should NOT have a _generate_descriptive_statistics_by_construct
# method.


# Function generalizes functionality from Roman's analyze_dual-crispr_NA_combined_simple-null-w-lfdr.R line 27 & 29,
# shown below with additional comments by ab:
#    # ab: x1 is log2 frequencies for the 1st replicate of all timepts
#    # ab: ab1 is abundance thresholds for all 1st replicates
#    good1 < -t(t(x1) > ab1)
#    # ab: g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the
#    # ab: first replicate.
#    # ab: 2 = min_num_timepoints_above_abundance_threshold
#    useless1 < -apply(good1, 1, sum) < 2
# Note that the R code does not parameterize min_num_timepoints_above_abundance_threshold but this code does.
def _get_constructs_to_ignore_mask_series(log2_fractions_by_constructs_by_timepoints_df,   # i.e., x1
                                          log2_fractions_thresholds_by_timepoints_series,  # i.e., ab1
                                          min_num_timepoints_above_abundance_threshold):
    """

    Args:
        min_num_timepoints_above_abundance_threshold (integer):

    Returns:
        pandas.Series:
    """
    passes_threshold_bool_by_construct_by_timept_df = _get_passes_threshold_bool_by_construct_by_timept_df(
        log2_fractions_by_constructs_by_timepoints_df, log2_fractions_thresholds_by_timepoints_series)
    num_timepts_passing_threshold_by_construct_series = passes_threshold_bool_by_construct_by_timept_df.sum(
        axis=1)  # axis of 1 or ‘columns’ means (unintuitively!) apply function across all column in a *row*--i.e.,
    # add up all the values in the row.  sum adds up Trues
    fails_min_num_timepts_bool_by_construct_series = num_timepts_passing_threshold_by_construct_series < \
                                                     min_num_timepoints_above_abundance_threshold
    return fails_min_num_timepts_bool_by_construct_series


# Function generalizes functionality from Roman's analyze_dual-crispr_NA_combined_simple-null-w-lfdr.R line 27 & 31,
# shown below with additional comments by ab:
#    # ab: x1 is log2 frequencies for the 1st replicate of all timepts
#    # ab: ab1 is abundance thresholds for all 1st replicates
#    good1 < -t(t(x1) > ab1)
#    # ab: g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the
#    # ab: first replicate.
#    # ab: useless1 = fails_min_num_timepts_bool_by_construct_series
#    good1[useless1,]<-FALSE #remove singletons
# Note that the R code does not parameterize min_num_timepoints_above_abundance_threshold but this code does.
def _get_usable_constructs_by_timepoints_mask_df(log2_fractions_by_constructs_by_timepoints_df,
                                                 log2_fractions_thresholds_by_timepoints_series,
                                                 fails_min_num_timepts_bool_by_construct_series):
    # for constructs that didn't have gte the minimum number of timepoints above the threshold, set FALSE for
    # all timepoints for that construct in the passes_threshold_bool_by_construct_by_timept_df (even for timepoints
    # that themselves had values above the threshold)

    passes_threshold_bool_by_construct_by_timept_df = _get_passes_threshold_bool_by_construct_by_timept_df(
        log2_fractions_by_constructs_by_timepoints_df, log2_fractions_thresholds_by_timepoints_series)
    passes_threshold_bool_by_construct_by_timept_df[fails_min_num_timepts_bool_by_construct_series.values] = False
    return passes_threshold_bool_by_construct_by_timept_df


# Function generalizes functionality from Roman's analyze_dual-crispr_NA_combined_simple-null-w-lfdr.R line 27,
# shown below with additional comments by ab:
#    # ab: x1 is log2 frequencies for the 1st replicate of all timepts
#    # ab: ab1 is abundance thresholds for all 1st replicates
#    # ab: constructs as rows, timept for 1st replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2
#    # ab: freq for row/col combination is above relevant abundance threshold
#    good1<-t(t(x1)>ab1)
def _get_passes_threshold_bool_by_construct_by_timept_df(log2_fractions_by_constructs_by_timepoints_df,    # i.e., x1
                                                         log2_fractions_thresholds_by_timepoints_series):  # i.e, ab1
    """
    Get construct x timepoint dataframe of booleans of whether each combination was above abundance threshold

    Returns:
        pandas.DataFrame:
    """
    passes_threshold_bool_by_construct_by_timept_df = log2_fractions_by_constructs_by_timepoints_df.ge(
        log2_fractions_thresholds_by_timepoints_series, axis="columns")
    return passes_threshold_bool_by_construct_by_timept_df


# Function generalizes functionality from Roman's analyze_dual-crispr_NA_combined_simple-null-w-lfdr.R lines 52-57,
# shown below with additional comments by ab:
#       # ab: i = construct index
#       # ab: g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the
#       # ab: first replicate
#       # ab: x1 is log2 frequencies for the 1st replicate of all timepts
#       # ab: if there are at least two good timepoints for this construct in replicate 1
#       if (sum(g1)>1) { #it's a good experiment
#          # ab: get the mean of the log2 frequencies for all the good timepoints for this construct in replicate 1
#          mx1<-mean(x1[i,g1])
#          # ab: get mean of timepoints (e.g. mean number of days) for all the good timepoints for this construct in
#          # ab: replicate 1
#          mt1<-mean(time[g1]) # ab: time is GLOBAL vector of timepoints (numbers, usually in days, in ascending order)
#          # ab: get variance of timepoints (e.g. mean number of days) for all the good timepoints for this construct in
#          # ab: replicate 1
#          v1<-Var(time[g1])
#          # ab: f1 = covariance of log2 frequencies for all the good timepoints for this construct in replicate 1 with
#          # ab: the timepoints for those good timepoints
#          f1<-Cov(x1[i,g1],time[g1])
#       }
# Note that the code above performs operations for a single construct at a time, while the code here calculates these
# values for all constructs at the same time.
def _generate_descriptive_statistics_by_construct(log2_fractions_by_constructs_by_timepoints_df,   # i.e., x1
                                                  usable_constructs_by_timepoint_mask):            # i.e., g1
    """

    Returns:
        pandas.DataFrame:
    """
    usable_log2_fractions_by_construct_df = log2_fractions_by_constructs_by_timepoints_df[
        usable_constructs_by_timepoint_mask]
    mean_usable_log2_fraction_by_construct_series = usable_log2_fractions_by_construct_df.mean(axis=1)
    # mx1 < - mean(x1[i, g1])

    # get mean of timepoints (e.g. mean number of days) and variance of timepoints for all usable timepoints
    # for this construct in this replicate
    timepoints_df = pandas.DataFrame(index=log2_fractions_by_constructs_by_timepoints_df.index,
                                     columns=log2_fractions_by_constructs_by_timepoints_df.columns.values)
    for curr_timept_header in log2_fractions_by_constructs_by_timepoints_df.columns.values:
        curr_timept, _ = ns_prep.read_timepoint_and_replicate_from_standardized_count_header(curr_timept_header)
        timepoints_df[curr_timept_header] = curr_timept
    usable_timepoints_by_construct_df = timepoints_df[usable_constructs_by_timepoint_mask]
    mean_usable_timepoints_by_construct_series = usable_timepoints_by_construct_df.mean(axis=1)
    # ddof = 0 indicates we want the *population* variance, not the sample variance, so normalize by n not n-1
    variance_usable_timepoints_by_construct_series = usable_timepoints_by_construct_df.var(axis=1, ddof=0)
    # anything with variance of NaN now has variance of 0
    variance_usable_timepoints_by_construct_series.fillna(0, inplace=True)

    covariance_per_construct_series = _calculate_covariance_by_construct_series(
        usable_log2_fractions_by_construct_df, usable_timepoints_by_construct_df)

    series = [mean_usable_log2_fraction_by_construct_series, mean_usable_timepoints_by_construct_series,
              variance_usable_timepoints_by_construct_series, covariance_per_construct_series]
    descriptive_stats_df = pandas.concat(series, axis=1)
    descriptive_stats_df.columns = [
        ThresholdedPerReplicateData.get_mean_usable_log2_fraction_for_construct_header(),  # called mx# in Roman's code
        ThresholdedPerReplicateData.get_mean_of_usable_timepoints_header(),                # called mt# in Roman's code
        ThresholdedPerReplicateData.get_variance_of_usable_timepoints_header(),            # called v# in Roman's code
        ThresholdedPerReplicateData.get_covariance_of_log2_fractions_w_timepoints_header()]  # called f# in Roman's code
    return descriptive_stats_df


# Sadly, as of this writing, pandas.DataFrame.cov does NOT take a ddof parameter, meaning that there is no way to change
# its default behavior (which calculates the sample covariance, using n-1 as the normalization factor) to the desired
# behavior (calculating the population covariance, using n as the normalization factor).  I originally tried writing
# this functionality using apply() on a DataFrame, rather than explicitly looping over each construct, but the apply()
# version took FOUR times as long as this version :(
def _calculate_covariance_by_construct_series(usable_log2_fractions_by_construct_df, usable_timepoints_by_construct_df):
    cov_by_construct_series = pandas.Series(index=usable_log2_fractions_by_construct_df.index)
    # TODO: Add check that usable_log2_fractions_by_construct_df and usable_timepoints_by_construct_df have same index?
    for curr_probe_pair_id in usable_log2_fractions_by_construct_df.index.values:
        # NB: for some reason, numpy.vstack--which puts some number of arrays together one on top of the other to make
        # a matrix--requires all the arrays be put into a tuple first rather than using *args.
        # The outcome of vstack is a 2-by-num_timepoints matrix where first row is the good log2 fractions values
        # for the current construct and second row is the corresponding good timepoint values, e.g.:
        #               MYEXP_T0_A    MYEXP_T14_A    MYEXP_T28_A
        # log2 fraction -11.01        -14.29         -15.99
        # timepoint     0             14             28
        matrix_of_features_rows_by_observations_cols = numpy.vstack(
            (usable_log2_fractions_by_construct_df.loc[curr_probe_pair_id],
             usable_timepoints_by_construct_df.loc[curr_probe_pair_id]))

        # covariance of log2 fractions for all good timepoints for this construct in this replicate
        # with the timepoint for those good timepoints.  Note the importance of setting the "bias" parameter to true,
        # because we want to normalize by the number of observations (because we want the *population* covariance), not
        # by the default number-of-observations-minus-one (which gives *sample* covariance)
        covariance_matrix = numpy.cov(matrix_of_features_rows_by_observations_cols, bias=True)

        # a "covariance matrix" actually has the variances of each variable down the diagonals, and the covariance of
        # each pair of variables in the off-diagonal cells--and it is symmetrical, because cov(x,y) = cov(y,x).  We
        # will only ever have two variables--log2 fractions and timepoints--so there will only be two off-diagonal
        # values, and they will always be the same (so take whichever one you want :)
        cov_by_construct_series[curr_probe_pair_id] = covariance_matrix[0, 1]

    cov_by_construct_series.fillna(0, inplace=True)  # set anything with covariance of NaN to have covariance of 0
    return cov_by_construct_series

# TODO: Decide if I can get rid of these implementations--no longer needed?
# def _calculate_descriptive_statistics_for_construct_info(construct_info_series,
#                                                          fails_min_num_timepts_bool_by_construct_series,
#                                                          usable_constructs_by_timepoint_mask):
#
#     # if there are at least min_num_timepoints_above_abundance_threshold good timepoints for this construct in
#     # this replicate:
#     if fails_min_num_timepts_bool_by_construct_series.loc[construct_info_series.name]:
#         stats_tuple = (None, None, 0, 0)
#     else:
#         stats_tuple = _calculate_descriptive_statistics_for_usable_construct_info(
#             construct_info_series, usable_constructs_by_timepoint_mask)
#
#     result = pandas.Series(data=stats_tuple,
#                            index=[ThresholdedPerReplicateData.get_mean_usable_log2_fraction_for_construct_header(),
#                                   ThresholdedPerReplicateData.get_mean_of_usable_timepoints_header(),
#                                   ThresholdedPerReplicateData.get_variance_of_usable_timepoints_header(),
#                                   ThresholdedPerReplicateData.get_covariance_of_log2_fractions_w_timepoints_header()],
#                            name=construct_info_series.name)
#     return result
#
#
# def _calculate_descriptive_statistics_per_construct(construct_info_series,
#                                                                 usable_constructs_by_timepoint_mask):
#
#     usable_timepoints_for_construct_mask = usable_constructs_by_timepoint_mask.loc[probe_pair_id]
#
#     # get the mean of the log2 fractions for all the good timepoints for this construct in this replicate
#     usable_log2_fractions_for_construct = construct_info_series[usable_timepoints_for_construct_mask]
#     mean_usable_log2_fraction_for_construct = numpy.mean(usable_log2_fractions_for_construct)
#     # mx1 < - mean(x1[i, g1])
#
#     # get mean of timepoints (e.g. mean number of days) and variance of timepoints for all usable timepoints
#     # for this construct in this replicate
#     usable_timepoints = ThresholdedPerReplicateData.extract_values_set_from_data_headers(
#         usable_log2_fractions_for_construct.to_frame().T, get_timepts=True)
#     mean_of_usable_timepoints = numpy.mean(usable_timepoints)
#     variance_of_usable_timepoints = numpy.var(usable_timepoints)
#
#     # NB: for some reason, numpy.vstack--which puts some number of arrays together one on top of the other to make
#     # a matrix--requires all the arrays be put into a tuple first rather than using *args.
#     # The outcome of vstack is a 2-by-num_timepoints matrix where first row is the good log2 fractions values
#     # for the current construct and second row is the corresponding good timepoint values, e.g.:
#     #               MYEXP_T0_A    MYEXP_T14_A    MYEXP_T28_A
#     # log2 fraction -11.01        -14.29         -15.99
#     # timepoint     0             14             28
#     matrix_of_features_rows_by_observations_cols = numpy.vstack(
#         (usable_log2_fractions_for_construct, usable_timepoints))
#
#     # covariance of log2 fractions for all good timepoints for this construct in this replicate
#     # with the timepoint for those good timepoints.  Note the importance of setting the "bias" parameter to true,
#     # because we want to normalize by the number of observations (because we want the *population* covariance), not
#     # by the default number-of-observations-minus-one (which gives *sample* covariance)
#     covariance_matrix = numpy.cov(matrix_of_features_rows_by_observations_cols, bias=True)
#
#     # a "covariance matrix" actually has the variances of each variable down the diagonals, and the covariance of
#     # each pair of variables in the off-diagonal cells--and it is symmetrical, because cov(x,y) = cov(y,x).  We
#     # will only ever have two variables--log2 fractions and timepoints--so there will only be two off-diagonal
#     # values, and they will always be the same (so take whichever one you want :)
#     covariance_of_log2_fractions_w_timepoints = covariance_matrix[0, 1]
#
#     result = (mean_usable_log2_fraction_for_construct, mean_of_usable_timepoints,
#               variance_of_usable_timepoints, covariance_of_log2_fractions_w_timepoints)
#     return result
