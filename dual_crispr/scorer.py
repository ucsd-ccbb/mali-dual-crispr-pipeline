# third-party libraries
import numpy
import pandas

# ccbb libraries
import dual_crispr.thresholded_per_replicate_data as ns_thresholded
import dual_crispr.fitted_per_replicate_data as ns_fitted

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


def temp_fit_ac_fc(per_replicate_data_list, min_num_timepoints_above_abundance_threshold):
    # get descriptive stats for each rep
    thresholded_per_rep_data_list = []
    for curr_per_replicate_data in per_replicate_data_list:
        curr_thresholded_per_rep_data = ns_thresholded.generate_thresholded_per_rep_data_from_per_rep_data(
            curr_per_replicate_data, min_num_timepoints_above_abundance_threshold)
        thresholded_per_rep_data_list.append(curr_thresholded_per_rep_data)

    naive_fitness_per_construct_series = _generate_naive_fitness_per_construct_series(thresholded_per_rep_data_list)

    fitted_per_rep_data_list = []
    for curr_thresholded_per_rep_data in thresholded_per_rep_data_list:
        curr_fitted_per_rep_data = ns_fitted.generate_fitted_per_rep_data_from_thresholded_per_rep_data(
            curr_thresholded_per_rep_data, naive_fitness_per_construct_series)

        # now do sdfc, tstat, pp_fc calculations


# def _get_constructs_to_ignore_across_replicates_mask(per_replicate_data_list,
#                                                      min_num_timepoints_above_abundance_threshold):
#     constructs_to_ignore_across_replicates = None
#     for curr_per_replicate_info in per_replicate_data_list:
#         curr_constructs_to_ignore_in_replicate = curr_per_replicate_info._get_constructs_to_ignore_mask_series(
#             min_num_timepoints_above_abundance_threshold)
#         if constructs_to_ignore_across_replicates is None:
#             constructs_to_ignore_across_replicates = curr_constructs_to_ignore_in_replicate
#         else:
#             constructs_to_ignore_across_replicates = constructs_to_ignore_across_replicates & \
#                                                      curr_constructs_to_ignore_in_replicate
#
#     return constructs_to_ignore_across_replicates


def _generate_naive_fitness_per_construct_series(thresholded_per_rep_data_list):
    descriptive_stats_per_construct_by_rep_num = _make_multiindex_df_across_replicates(
        thresholded_per_rep_data_list, "descriptive_stats_df")

    sum_covariances_per_construct_across_replicates_series = _calc_sum_of_stat_per_construct_across_replicates_series(
        descriptive_stats_per_construct_by_rep_num,
        ns_thresholded.ThresholdedPerReplicateData.get_covariance_of_log2_fractions_w_timepoints_header())
    sum_variances_per_construct_across_replicates_series = _calc_sum_of_stat_per_construct_across_replicates_series(
        descriptive_stats_per_construct_by_rep_num,
        ns_thresholded.ThresholdedPerReplicateData.get_variance_of_usable_timepoints_header())

    # fc[i] is the combined fitness (across replicates) for construct i
    # fc[i] < - (f1 + f2) / (v1 + v2)  # the combined fitness from replicate 1+2
    fitness_per_construct_across_replicates_series = sum_covariances_per_construct_across_replicates_series.div(
        sum_variances_per_construct_across_replicates_series)
    fitness_per_construct_across_replicates_series.fillna(0, inplace=True)  # NaNs will become 0 per Roman's code
    return fitness_per_construct_across_replicates_series


def _generate_stddev_of_naive_fitness_per_construct_series(fitted_per_rep_data_list):
    log2_fract_per_construct_per_timept_per_replicate_df = _make_multiindex_df_across_replicates(
        fitted_per_rep_data_list, "_log2_fractions_by_constructs_by_timepoints_df")
    expected_log2_fract_per_construct_per_timept_per_replicate_df = _make_multiindex_df_across_replicates(
        fitted_per_rep_data_list, "_expected_log2_fraction_per_construct_per_timepoint_df")
    usable_constructs_by_timepoint_mask_per_replicate = _make_multiindex_df_across_replicates(
        fitted_per_rep_data_list, "_usable_constructs_by_timepoint_mask")
    concatenated_timepoints_series_per_replicate = _make_multiindex_df_across_replicates(fitted_per_rep_data_list,
                                                                                          "timepoints_series", axis=0)

    # sdfc = numpy array len of num of constructs filled with 0.1
    num_constructs = len(log2_fract_per_construct_per_timept_per_replicate_df.index.values)
    sdfc = numpy.array([1, num_constructs])
    sdfc.fill(0.1)

    # for each construct
    for curr_index in log2_fract_per_construct_per_timept_per_replicate_df.index.values:
        pass
        #stat_per_construct_across_replicates_df = stats_per_construct_by_rep_num_df.xs(stat_header, level=1, axis=1)
        # mask = _usable_constructs_by_timepoint_mask_per_replicate[currConstruct]
        # if there are no TRUEs anywhere in mask
            # currxfits = xfits[currConstruct, mask]
            # curr_log2_fractions = log2 fractions [currConstruct, mask]
            # numerator_dif = currxfits - curr_log2_fractions
            # numerator = sqrtsum(numerator_dif)

            # curr_timepts = concatenated_timepts_series[mask]
            # mean_curr_timepts = mean(curr_timepts)
            # denominator_diff = curr_timepts - mean_curr_timepts
            # denominator = sqrtsum(denominator_dif)

            # sdfc[currConstruct] = numerator / denominator
        # end if
    # next construct


#
# # per construct:
#     # get expected fits for all the good replicates for this construct
#     # get log2 frequencies for all the good replicates for this construct
#     # subtract log2 frequencies vector from expected fits vector
#     # calc sqrtsum of difference vector
#
#     # get the timepoints for all the good replicates for this construct
#     # find the mean of the above
#     # subtract the mean from each entry in the timepoints vector
#     # call sqrtsum of difference vector
#
#     # divide first difference vector by second difference vector
#
#     sdfc < - rep(0.1, nx)  # standard error of fc
#
#     # for 1 to number of constructs
#     for (i in 1: nx) {
#     # if this construct doesn't have enough "good" measurements, skip it
#     if (allbad[i])
#     next
#
#     # g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the first replicate
#     g1 < - good1[i,]
#     g2 < - good2[i,]
#
#     # sqrtsum<-function(y) sqrt(sum(y^2))
#     # outside of this function, this sdfc value is used only in the output of the construct file
#     # I think that sdfc is "standard deviation of fc" for each construct.
#
#     # So, Roman says:
#     # std err of fc = sqrt of sum over t of (lower-case-epsilon for construct c, as a function of t)^2
#     # divided by sqrt of (nc - 2)*sum over t of (t^2 - (mean of t)^2)
#     # Note that lower-case-epsilon is Xc(t) - xc(t) = xX - xfitX
#     # so xfitX - xX = - lower-case-epsilon ... but since it is being squared and then square-rooted, I
#     # suppose the negation doesn't matter.
#     # nc - 2 is the number of degrees of freedom
#     # where nc = number of data points = 2*nt minus any number of points below the threshold
#     # (note the description above seems to assume 2 replicates) ... seems to me this must mean
#     # number of data points *for this construct c*, not total.
#     # nt in above is number of timepoints, as here ...
#     # Roman also says that tc [i.e., the t statistic for construct c] = fc /SE(fc)
#     # and the internet tells me that SE(x) = SD(x)/sqrt(n) ... but maybe the sqrt(n) is just the most usual
#     # sqrt of degrees of freedom, and could be something else in a more complex system.
#     # So, the value being calculated directly below is SD(x), which is why it doesn't have the
#     # nc - 2 term that is in the denominator of the SE(x) calculation (in the manuscript, Roman says that
#     # fc's sd = sqrt(nc-2)*SE(fc), where nc-2 = degrees of freedom.  Farther below, after the calculation
#     # of sdfc, we get the calculation of tc, which is fc/(sdfc/sqrt(df)), where the sdfc/sqrt(df) term
#     # is the calculation of SE(fc).
#
#     # result is vector with one std dev of fc for each construct
#     sdfc[i] < -
#               sqrtsum(c(xfit1[i, g1], xfit2[i, g2]) - c(x1[i, g1], x2[i, g2])) /
#               sqrtsum(c(time[g1], time[g2]) - mean(c(time[g1], time[g2])))
#
#     }


def _make_multiindex_df_across_replicates(per_rep_obj_list, name_of_df_attribute, axis=1):
    """

    Args:
        thresholded_per_rep_data_list (list):

    Returns:
        pandas.DataFrame
    """
    level_one_index_names = [x.replicate_num for x in per_rep_obj_list]
    list_of_dfs = [getattr(x, name_of_df_attribute) for x in per_rep_obj_list]
    multiindex_df = pandas.concat(list_of_dfs, axis=axis, keys=level_one_index_names)
    return multiindex_df


def _calc_sum_of_stat_per_construct_across_replicates_series(stats_per_construct_by_rep_num_df, stat_header):
    stat_per_construct_across_replicates_df = stats_per_construct_by_rep_num_df.xs(stat_header, level=1, axis=1)
    sum_stat_per_construct_across_replicates_series = stat_per_construct_across_replicates_df.sum(axis=1)
    return sum_stat_per_construct_across_replicates_series
