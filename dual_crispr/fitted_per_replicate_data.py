# third-party libraries
import numpy
import pandas

# ccbb libraries
from dual_crispr.thresholded_per_replicate_data import ThresholdedPerReplicateData

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class FittedPerReplicateData(ThresholdedPerReplicateData):
    # NB: _run_init_methods param only exists to facilitate unit-testing and should NOT be changed for any other usage
    def __init__(self, replicate_num, log2_fractions_by_constructs_by_timepoints_df,
                 log2_fractions_thresholds_by_timepoints_series, min_num_timepoints_above_abundance_threshold,
                 unnormed_fitness_per_construct_series, _run_init_methods=True):

        super().__init__(replicate_num, log2_fractions_by_constructs_by_timepoints_df,
                         log2_fractions_thresholds_by_timepoints_series, min_num_timepoints_above_abundance_threshold)

        if _run_init_methods:
            self._normed_init_log2_fraction_per_construct_series = \
                _calculate_normed_init_log2_fraction_per_construct_series(
                    self._descriptive_stats_by_construct_df, unnormed_fitness_per_construct_series,
                    self._log2_fractions_by_constructs_by_timepoints_df,
                    self._fails_min_num_timepts_bool_by_construct_series)
            self._lambda_per_timepoint_series = _calculate_lambda_per_timepoint_series(
                self._normed_init_log2_fraction_per_construct_series,unnormed_fitness_per_construct_series,
                self.timepoints_series)
            self._expected_log2_fraction_per_construct_per_timepoint_df = \
                _calculate_expected_log2_fraction_per_construct_per_timepoint_df(
                    self.timepoints_series,
                    self._normed_init_log2_fraction_per_construct_series, unnormed_fitness_per_construct_series,
                    self._lambda_per_timepoint_series)


# Not a class method of FittedPerReplicateData because anything that inherits from FittedPerReplicate Data shouldn't
# have this method
def generate_fitted_per_rep_data_from_thresholded_per_rep_data(
        thresholded_per_replicate_data, unnormed_fitness_per_construct_series):
    """

    Args:
        unnormed_fitness_per_construct_series (pandas.Series):
        thresholded_per_replicate_data (ThresholdedPerReplicateData):
    """
    result = FittedPerReplicateData(thresholded_per_replicate_data.replicate_num,
                                    thresholded_per_replicate_data._log2_fractions_by_constructs_by_timepoints_df,
                                    thresholded_per_replicate_data._log2_fractions_thresholds_by_timepoints_series,
                                    thresholded_per_replicate_data._min_num_timepoints_above_abundance_threshold,
                                    unnormed_fitness_per_construct_series)
    return result


# Note: the functions below are NOT methods of the FittedPerReplicateData object--even though that is the only
# entity that uses their functionality--because they should not be available to a FittedPerReplicateData object
# after it is created.  That is, once you've run _calculate_expected_log2_fraction_per_construct_per_timepoint_df,
# you should not be running it again on the FittedPerReplicateData object later.


# Function generalizes functionality from Roman's MethodII.R lines 38, 68-70, and 75-76, shown below with additional
# comments by ab:
#    ac1<-x1[,1] #just a guess # ab: log2 frequencies for all constructs for first timepoint in this replicate
#
#       # ab: mx1 = mean of the log2 frequencies for all the good timepoints for this construct in replicate 1
#       # ab: mt1 = mean of timepoints (e.g. mean number of days) for all the good timepoints for this construct in
#       # ab: replicate 1
#       # ab: fc[i] is the combined fitness (across replicates) for construct i
#       if (sum(g1)>1) {
#           # ab: if there are at least two good timepoints for this construct in replicate 1
#           # ab: ac1 = mean of log2 freqs for this construct for good timepoints for rep 1 - (mean of timepts for good
#           # ab: timepoints for this construct for rep 1)*combined fitness across replicates for this construct
#          ac1[i]<-mx1-fc[i]*mt1
#       }
#
#    # ab: ac is the initial condition (in log2 frequency) for each construct
#    # ab: this is normalizing ac1:
#    # ab: Roman's methods document say "By definition, log2 relative frequencies satisfy the constraint
#    # ab: sum over c of (2^xc) = 1 at all times"
#    # ab: ac1 is the set of log2 frequencies for all constructs in replicate 1 at "initial conditions", so
#    # ab: ac1 must satisfy the above constraint.  IFF the above constraint is satisfied, -log2(1) = 0.
#    # ab: If the above constraint isn't satisfied, the value of -log2(sum over c of (2^ac)) is subtracted from
#    # ab: ac to ensure the constraint is satisfied.
#    alpha<- -log2(sum(2^ac1))
#    ac1<-ac1+alpha #enforce normalization at time=0, sum(2^ac)=1
def _calculate_normed_init_log2_fraction_per_construct_series(descriptive_stats_per_construct,
                                                             unnormed_fitness_per_construct_series,
                                                             log2_fractions_by_constructs_by_timepoints_df,
                                                             fails_min_num_timepts_bool_by_construct_series):

    mean_log2_header = ThresholdedPerReplicateData.get_mean_usable_log2_fraction_for_construct_header()
    mean_timept_header = ThresholdedPerReplicateData.get_mean_of_usable_timepoints_header()

    # ac is the initial condition (in log2 frequency) for construct c
    # ac1 of each construct = mean of log2 freqs for this construct for good timepoints for rep 1
    # - [(mean of timepts for good timepoints for each construct for rep 1)
    #    * combined fitness across replicates for each construct]
    unnormed_init_log2_fraction_per_construct_series = descriptive_stats_per_construct[mean_log2_header] \
        - (unnormed_fitness_per_construct_series * descriptive_stats_per_construct[mean_timept_header])

    # for constructs lacking at least min_num_timepoints_above_abundance_threshold good timepoints in this replicate:
    # just guess: set starting value of init_log2_fraction to log2 frequencies for construct for first timepoint
    # in this replicate
    failing_constructs_mask = fails_min_num_timepts_bool_by_construct_series.values
    unnormed_init_log2_fraction_per_construct_series[failing_constructs_mask] = \
        log2_fractions_by_constructs_by_timepoints_df.ix[failing_constructs_mask, 0]

    # Now normalize ac1:
    # Roman's methods say "By definition, log2 relative frequencies satisfy the constraint
    # sum over c of (2^xc) = 1 at all times"
    unnormed_init_fraction_per_construct_series = numpy.exp2(unnormed_init_log2_fraction_per_construct_series)
    norm_const = -numpy.log2(unnormed_init_fraction_per_construct_series.sum())
    normed_init_log2_fraction_per_construct_series = unnormed_init_log2_fraction_per_construct_series + norm_const
    return normed_init_log2_fraction_per_construct_series  # i.e., ac1


# Function generalizes functionality from Roman's MethodII.R lines 80-83, shown below with additional comments by ab:
#    # ab: ac1 is the initial condition (in log2 frequency) for each construct in replicate 1
#    # ab: fc is the combined fitness (across replicates) for each construct
#    # ab: time is GLOBAL vector of timepoints (numbers, usually in days, in ascending order)
#    # ab: for 1 to number of timepoints
#    for (i in 1:nt) {
#       lambda1[i]<- -log2(sum(2^(ac1+fc*time[i])))
#       ...
#    } #these are initial estimates of lambda(t)
# Note that the code above performs operations for a single timepoint at a time, while the code here calculates these
# values for all timepoints at the same time.
def _calculate_lambda_per_timepoint_series(normed_init_log2_fraction_per_construct_series,             # i.e., ac1
                                           unnormed_fitness_per_construct_series, timepoints_series):  # i.e., fc, time

    # NB: timepoints are NOT limited to "usable" timepoints, per Roman's original implementation
    exponent_per_construct_per_timepoint_df = \
        _calculate_normed_init_log2_fraction_plus_unnormed_fitness_mult_by_time_per_construct_per_timept_df(
            normed_init_log2_fraction_per_construct_series, unnormed_fitness_per_construct_series, timepoints_series)
    two_to_exponent_per_construct_per_timepoint_df = numpy.exp2(exponent_per_construct_per_timepoint_df)
    sum_across_constructs_per_timepoint_series = two_to_exponent_per_construct_per_timepoint_df.sum(axis=0)
    lambda_per_timepoint_series = -numpy.log2(sum_across_constructs_per_timepoint_series)
    return lambda_per_timepoint_series  # i.e., lambda1


# Function generalizes calculation ac1+fc*time[i] from Roman's MethodII.R lines 81 & 87.
#    # ab: i = current *timepoint* (not construct)
#    # ab: ac1 is the initial condition (in log2 frequency) for each construct in replicate 1
#    # ab: fc is the combined fitness (across replicates) for each construct
#    # ab: time is GLOBAL vector of timepoints (numbers, usually in days, in ascending order)
# Note that the code above performs operations for a single timepoint at a time, while the code here calculates these
# values for all timepoints at the same time.
def _calculate_normed_init_log2_fraction_plus_unnormed_fitness_mult_by_time_per_construct_per_timept_df(
        normed_init_log2_fraction_per_construct_series, unnormed_fitness_per_construct_series, timepoints_series):
    # i.e., inputs are ac1, fc, and time

    unnormed_fitness_mult_by_time_per_construct_per_timept_df = \
        _calculate_unnormed_fitness_mult_by_time_per_construct_per_timept_df(
            unnormed_fitness_per_construct_series, timepoints_series)

    normed_init_log2_fraction_plus_unnormed_fitness_mult_by_time_per_construct_per_timept_df = \
        unnormed_fitness_mult_by_time_per_construct_per_timept_df.add(normed_init_log2_fraction_per_construct_series,
                                                                   axis='index')
    return normed_init_log2_fraction_plus_unnormed_fitness_mult_by_time_per_construct_per_timept_df  # i.e., ac1+fc*time


# Function generalizes calculation fc*time[i] from Roman's MethodII.R lines 81 & 87.
#    # ab: i = current *timepoint* (not construct)
#    # ab: fc is the combined fitness (across replicates) for each construct
#    # ab: time is GLOBAL vector of timepoints (numbers, usually in days, in ascending order)
# Note that the code above performs operations for a single timepoint at a time, while the code here calculates these
# values for all timepoints at the same time.
def _calculate_unnormed_fitness_mult_by_time_per_construct_per_timept_df(unnormed_fitness_per_construct_series,
                                                                      timepoints_series):  # i.e, fc and time
    # NB: per Roman's original implementation, timepoints are NOT limited to "usable" timepoints
    # Yes, it would be simpler to use series.mul, but also >1 order of magnitude slower
    numpy_matrix = unnormed_fitness_per_construct_series.values[:, numpy.newaxis] * timepoints_series.values
    unnormed_fitness_mult_by_time_per_construct_per_timept_df = pandas.DataFrame(
        numpy_matrix, columns=timepoints_series.index.values, index=unnormed_fitness_per_construct_series.index)
    return unnormed_fitness_mult_by_time_per_construct_per_timept_df  # i.e., fc*time


# Function generalizes functionality from Roman's MethodII.R lines 85-88, shown below with additional comments by ab:
#    # ab: x1 is log2 frequencies for the 1st replicate of all timepts
#    xfit1<-x1 #for size
#    # ab: I think the "for size" comment means that xfitX is being set to xX not because we're using *any* of the xX
#    # ab: values but just to initialize xfitX to the desired size (which is the same as the size of xX)
#    # ab: ac1 is the initial condition (in log2 frequency) for each construct in replicate 1
#    # ab: fc is the combined fitness (across replicates) for each construct
#    # ab: time is GLOBAL vector of timepoints (numbers, usually in days, in ascending order)
#    # ab: nt = number of timepoints
#    # ab: lambda1[i] = -log2(sum(2^(ac1+fc*time[i])))
#    for (j in 1:nt) {
#       # ab: the expected value of xc(t) at this timepoint t, calculated as a column for all constructs c
#       xfit1[,j]<-ac1+fc*time[j]+lambda1[j]
#    }
# Note that the code above performs operations for a single timepoint at a time, while the code here calculates these
# values for all timepoints at the same time.
def _calculate_expected_log2_fraction_per_construct_per_timepoint_df(
        timepoints_series, normed_init_log2_fraction_per_construct_series, unnormed_fitness_per_construct_series,
        lambda_per_timepoint_series):
    # i.e., inputs are time, ac1, fc, lambda1

    normed_init_log2_fraction_plus_unnormed_fitness_mult_by_time_per_construct_per_timept_df = \
        _calculate_normed_init_log2_fraction_plus_unnormed_fitness_mult_by_time_per_construct_per_timept_df(
            normed_init_log2_fraction_per_construct_series, unnormed_fitness_per_construct_series, timepoints_series)

    lambdas_transposed = lambda_per_timepoint_series.values.transpose()
    numpy_matrix = lambdas_transposed + \
                   normed_init_log2_fraction_plus_unnormed_fitness_mult_by_time_per_construct_per_timept_df.values
    expected_log2_fraction_per_construct_per_timepoint_df = pandas.DataFrame(
        numpy_matrix,
        columns=normed_init_log2_fraction_plus_unnormed_fitness_mult_by_time_per_construct_per_timept_df.columns.values,
        index=normed_init_log2_fraction_plus_unnormed_fitness_mult_by_time_per_construct_per_timept_df.index)

    return expected_log2_fraction_per_construct_per_timepoint_df  # i.e., xfit1
