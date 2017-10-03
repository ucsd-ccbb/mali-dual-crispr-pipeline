# standard libraries
import contextlib
import math
import os
import subprocess
import tempfile

# third-party libraries
import numpy
import pandas
import scipy.stats

# ccbb libraries
import dual_crispr.thresholded_per_replicate_data as ns_thresholded
import dual_crispr.fitted_per_replicate_data as ns_fitted

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


def get_num_constructs(any_per_rep_data_list):
    # TODO: add some sort of checking that all per_rep_data has same num of constructs
    # TODO: refactor reference to private instance property
    return len(any_per_rep_data_list[0]._log2_fractions_by_constructs_by_timepoints_df.index.values)


def get_num_timepoints(any_per_rep_data_list):
    # TODO: add some sort of checking that all per_rep_data has same num of timepoints
    return len(any_per_rep_data_list[0].timepoints_series)


# Function generalizes functionality from Roman's analyze_dual-crispr_NA_combined_simple-null-w-lfdr.R lines 22-130,
# shown below with additional comments by ab:
#
# ab: x1 is log2 frequencies for the 1st replicate of all timepts
# ab: x2 is log2 frequencies for the 2nd replicate of all timepts
# ab: ab1 is abundance thresholds for all 1st replicates
# ab: ab2 is abundance thresholds for all 2nd replicates
#
# fit_ac_fc<-function(x1,ab1,x2,ab2) { #badx is TRUE when x-value is bad
#
#    er_ac<-1
#    l<-0
#    nx<-nrow(x1)
#
# ----SEE IMPLEMENTATION IN: thresholded_per_replicate_data.py
#    # ab: constructs as rows, timept for 1st replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2
#    # ab: freq for row/col combination is above relevant abundance threshold
#    good1<-t(t(x1)>ab1)
#    # ab: constructs as rows, timept for 2nd replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2
#    # ab: freq for row/col combination is above relevant abundance threshold
#    good2<-t(t(x2)>ab2)
#    # ab: 1 = sum over rows--i.e., constructs.  Any construct that isn't above abundance threshold in 1st replicate in
#    # ab: at least 2 timepoints has 1/TRUE in "useless" matrix
#    useless1<-apply(good1,1,sum)<2
#    useless2<-apply(good2,1,sum)<2
#    # ab: For each constructs that isn't above abundance thresholds in at least two timepoints for this replicate, this
#    # ab: sets that construct's "good" values to FALSE for *all* timepoints in this replicate
#    good1[useless1,]<-FALSE #remove singletons
#    good2[useless2,]<-FALSE #remove singletons
# ----END SEE IMPLEMENTATION IN: thresholded_per_replicate_data.py
#
#    # ab: Note: here apply(goodX,1,sum) is NOT uselessX because goodX was changed above
#    # ab: allbad is true for all constructs that lack at least 2 acceptable-abundance timepoints in BOTH experiments
# #  allbad<-!(apply(good1[,1:2],1,sum)==2 | apply(good1[,(nt-1):nt],1,sum)==2) | !(apply(good2[,1:2],1,sum)==2 | apply(good2[,(nt-1):nt],1,sum)==2)
#    allbad<-apply(good1,1,sum)<2 & apply(good2,1,sum)<2 #in this case I have nothing to use in either experiment
#
#    # ab: nt = number of timepoints
#    lambda1<-rep(0,nt)
#    lambda2<-rep(0,nt)
#    ac1<-x1[,1] #just a guess # ab: log2 frequencies for all constructs for first timepoint in this replicate
#    ac2<-x2[,1] #just a guess
# #  ac1<-rep(-Inf,nx) #underrepresented constructs will get zero abundance
# #  ac2<-rep(-Inf,nx) #underrepresented constructs will get zero abundance
#    fc<-rep(0,nx)
#
# ----SEE IMPLEMENTATION IN: thresholded_per_replicate_data.py
#    # ab: for 1 to number of constructs
#    for (i in 1:nx) {
#       # ab: if this construct doesn't have at least two timepoints above abundance threshold in at least one
#       # ab: replicate, ignore it and move on
#       if (allbad[i]) next #from now on there is at least one good experiment
#
#       # ab: apparently "f" stands for fitness--which is covariance (see below)--and # stands for replicate
#       # ab: v stands for variance and # " "
#       f1<-0
#       f2<-0
#       v1<-0
#       v2<-0
#       # ab: g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the
#       # ab: first replicate
#       g1<-good1[i,]
#       # ab: if there are at least two good timepoints for this construct in replicate 1
#       if (sum(g1)>1) { #it's a good experiment
#          # ab: get the mean of the log2 frequencies for all the good timepoints for this construct in replicate 1
#          mx1<-mean(x1[i,g1])
#          # ab: get mean of timepoints (e.g. mean number of days) for all the good timepoints for this construct in
#          # ab: replicate 1
#          mt1<-mean(time[g1]) # ab: time is GLOBAL vector of timepoints (numbers, usually in days, in ascending order)
#          v1<-Var(time[g1]) # ab: get variance of timepoints " "
#          # ab: f1 = covariance of log2 frequencies for all the good timepoints for this construct in replicate 1 with
#          # ab: the timepoints for those good timepoints
#          f1<-Cov(x1[i,g1],time[g1])
#       }
#
#       # ab: do the exact same thing as above, but for replicate 2
#       g2<-good2[i,]
#       if (sum(g2)>1) { #it's a good experiment
#          mx2<-mean(x2[i,g2])
#          mt2<-mean(time[g2])
#          v2<-Var(time[g2])
#          f2<-Cov(x2[i,g2],time[g2])
#       }
# ----END SEE IMPLEMENTATION IN: thresholded_per_replicate_data.py
#
#       # ab: fc[i] is the combined fitness (across replicates) for construct i
#       fc[i]<-(f1+f2)/(v1+v2) #the combined fitness from replicate 1+2
#       #fc remains defined up to an additive constant
#
# ----SEE IMPLEMENTATION IN: fitted_per_replicate_data.py
#       if (sum(g1)>1) {
#           # ab: if there are at least two good timepoints for this construct in replicate 1
#           # ab: ac1 = mean of log2 freqs for this construct for good timepoints for rep 1 - (mean of timepts for good
#           # ab: timepoints for this construct for rep 1)*combined fitness across replicates for this construct
#          ac1[i]<-mx1-fc[i]*mt1
#       }
#       # ab: same as above but for replicate 2
#       if (sum(g2)>1) {
#          ac2[i]<-mx2-fc[i]*mt2
#       }
#    }
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
#    # ab: same as above but for replicate 2
#    alpha<- -log2(sum(2^ac2))
#    ac2<-ac2+alpha #enforce normalization at time=0, sum(2^ac)=1
#
#    # ab: ok, the *expected* log2 frequency for construct x at time t is:
#    # ab: xc(t) = ac + fc*t - log2(sum over c of 2^(ac + fc*t))
#    # ab: The below code calls the term starting with "-log2" "lambda".
#    # ab: Note that lambda is being calculated separately for each combination of timepoint+replicate
#
#    # ab: for 1 to number of timepoints
#    for (i in 1:nt) {
#       lambda1[i]<- -log2(sum(2^(ac1+fc*time[i])))
#       lambda2[i]<- -log2(sum(2^(ac2+fc*time[i])))
#    } #these are initial estimates of lambda(t)
#
#    # ab: x1 is log2 frequencies for the 1st replicate of all timepts
#    xfit1<-x1 #for size
#    # ab: I think the "for size" comment means that xfitX is being set to xX not because we're using *any* of the xX
#    # ab: values but just to initialize xfitX to the desired size (which is the same as the size of xX)
#    for (j in 1:nt) {
#       # ab: the expected value of xc(t) at this timepoint t, calculated as a column for all constructs c
#       xfit1[,j]<-ac1+fc*time[j]+lambda1[j]
#    }
#    # ab: same as above but for replicate 2
#    xfit2<-x2 #for size
#    for (j in 1:nt) {
#       xfit2[,j]<-ac2+fc*time[j]+lambda2[j]
#    }
# ----END SEE IMPLEMENTATION IN: fitted_per_replicate_data.py
#
#    sdfc<-rep(0.1,nx) #standard error of fc
#    tstat<-rep(0,nx) # ab: presumably t statistic
#    # ab: df = vector w one entry for each construct, with value from 0 to -2; 0 if both replicates of relevant
#    # ab: construct are good, -1 if just one is, -2 if neither are.  Inits to 0 for all
#    df<-rep(0,nx)
#    # ab: p value from t test ... initialize to 1 (not significant) for everything
#    p_t<-rep(1,nx)
# #  sefc<-apply(xfit-x,1,sqrtsum)/sqrtsum(time)
#    # ab: for 1 to number of constructs
#    for (i in 1:nx) {
#       # ab: if this construct doesn't have enough "good" measurements, skip it
#       if (allbad[i]) next
#
#       # ab: g1 = true/false values of whether construct i passes various abundance filters for all
#       # ab: timepoints for the first replicate; g2 is analogous
#       g1 <- good1[i,]
#       g2 <- good2[i,]
#       # ab: df = vector w one entry for each construct, with value from 0 to -2; 0 if both replicates of relevant
#       # ab: construct are good, -1 if just one is, -2 if neither are
#       # ab: I suspect that "df" is "degrees of freedom" for each construct
#       df[i]<-sum(g1)+sum(g2)-2
#
#       # ab: reminder: sqrtsum<-function(y) sqrt(sum(y^2))
#       # ab: outside of this function, this sdfc value is used only in the output of the construct file
#       # ab: I think that sdfc is "standard deviation of fc" for each construct.
#
#       # ab: So, Roman says:
#       # ab: std err of fc = sqrt of sum over t of (lower-case-epsilon for construct c, as a function of t)^2
#       # ab: divided by sqrt of (nc - 2)*sum over t of (t^2 - (mean of t)^2)
#       # ab: Note that lower-case-epsilon is Xc(t) - xc(t) = xX - xfitX
#       # ab: so xfitX - xX = - lower-case-epsilon ... but since it is being squared and then square-rooted, I
#       # ab: suppose the negation doesn't matter.
#       # ab: nc - 2 is the number of degrees of freedom
#       # ab: where nc = number of data points = 2*nt minus any number of points below the threshold
#       # ab: (note the description above seems to assume 2 replicates) ... seems to me this must mean
#       # ab: number of data points *for this construct c*, not total.
#       # ab: nt in above is number of timepoints, as here ...
#       # ab: Roman also says that tc [i.e., the t statistic for construct c] = fc /SE(fc)
#       # ab: and the internet tells me that SE(x) = SD(x)/sqrt(n) ... but maybe the sqrt(n) is just the most usual
#       # ab: sqrt of degrees of freedom, and could be something else in a more complex system.
#       # ab: So, the value being calculated directly below is SD(x), which is why it doesn't have the
#       # ab: nc - 2 term that is in the denominator of the SE(x) calculation (in the manuscript, Roman says that
#       # ab: fc's sd = sqrt(nc-2)*SE(fc), where nc-2 = degrees of freedom.  Farther below, after the calculation
#       # ab: of sdfc, we get the calculation of tc, which is fc/(sdfc/sqrt(df)), where the sdfc/sqrt(df) term
#       # ab: is the calculation of SE(fc).
#
#       # ab: result is vector with one std dev of fc for each construct
# #      sefc[i]<-sqrtsum( c(xfit1[i,good1[i,]],xfit2[i,good2[i,]]) - c(x1[i,good1[i,]],x2[i,good2[i,]]) ) /sqrtsum( c(time[good1[i,]],time[good2[i,]]) - mean(c(time[good1[i,]],time[good2[i,]])) ) /sqrt(df[i])
#       sdfc[i]<-sqrtsum( c(xfit1[i,g1],xfit2[i,g2]) - c(x1[i,g1],x2[i,g2]) ) /sqrtsum( c(time[g1],time[g2]) - mean(c(time[g1],time[g2])) )
#    }
# #what is median sd?
#    has_sd<-df>0
#    median_sd<-median(sdfc[has_sd])
#    sdfc[!has_sd]<-median_sd #just so it isn't 0
#
#
#    # ab: for 1 to number of constructs
#    for (i in 1:nx) {
#       # ab: don't try to calculate t statistic and p value for any fc that doesn't have a stderr
#       if (!has_sd[i]) next
#       # ab: calc t statistic of fc
#       tstat[i]<-fc[i]/(sdfc[i]/sqrt(df[i]))
#       # ab: pt is R function, presumably to do t-test :)  Suppose pt result must be one-tailed, hence * 2
#       p_t[i]<-2*pt(-abs(tstat[i]),df=df[i]) #raw p-values from t-test
#    }
#
#    # ab: ok, lfdr stands for "local fdr", and the lfdr function comes from the qvalue bioconductor package
#    # ab: that is installed way above.  The first input is the vector of p-values, and the second input is the
#    # ab: estimated proportion of true null p-values, where pi0.method is "the method for automatically choosing tuning
#    # ab: parameter in the estimation of π0, the proportion of true null hypotheses."
#    # ab: nx = number of constructs; default value of lfdr is set to one for all of them
#    lfdr_fc<-rep(1,nx)
#    # ab: for constructs that have an sd, the lfdr of the fc could be calculated and is now set to its calculated value
#    # ab:  instead of the default
#    l<-lfdr(p_t[has_sd],pi0.method="bootstrap")
#    lfdr_fc[has_sd]<-l
#
#
# #  ptraw<-2*pt(-abs(tstat),df=df)
# #  lfdrt<-lfdr(ptraw,pi0.method="bootstrap")
#
#    # ab: ac1 is the initial condition (in log2 frequency) for each construct c for replicate 1
#    # ab: ac2 is the initial condition (in log2 frequency) for each construct c for replicate 2
#    # ab: fc is the fitness of each construct c (calculated across both replicates)
#    # ab: sdfc is the std deviation of the fitness of each construct c (calculated across both replicates)
#    # ab: p_t is the raw p value of the fc of each construct c (calculated across both replicates)
#    # ab: lfdr_fc is the local FDR of each construct (calculated across both replicates)
#    # ab: df is the degrees of freedom of each construct c (calculated across both replicates)
#    # ab: allbad is a boolean value for each construct c that is true for all the constructs that lack at least 2
#    # ab: acceptable-abundance timepoints in BOTH experiments
#    vl<-list(ac1,ac2,fc,sdfc,p_t,lfdr_fc,df,allbad)
#    return(vl)
# }
# TODO: Rename function to better reflect what it does ...
def temp_fit_ac_fc(per_replicate_data_list, min_num_timepoints_above_abundance_threshold):
    # get descriptive stats for each rep
    thresholded_per_rep_data_list = []
    for curr_per_replicate_data in per_replicate_data_list:
        curr_thresholded_per_rep_data = ns_thresholded.generate_thresholded_per_rep_data_from_per_rep_data(
            curr_per_replicate_data, min_num_timepoints_above_abundance_threshold)
        thresholded_per_rep_data_list.append(curr_thresholded_per_rep_data)

    constructs_to_ignore_across_replicates_series = _get_constructs_to_ignore_across_replicates_mask(
        thresholded_per_rep_data_list)

    unnormed_fitness_per_construct_series = _generate_unnormed_fitness_per_construct_series(
        thresholded_per_rep_data_list)

    fitted_per_rep_data_list = []
    for curr_thresholded_per_rep_data in thresholded_per_rep_data_list:
        curr_fitted_per_rep_data = ns_fitted.generate_fitted_per_rep_data_from_thresholded_per_rep_data(
            curr_thresholded_per_rep_data, unnormed_fitness_per_construct_series)
        fitted_per_rep_data_list.append(curr_fitted_per_rep_data)

    stddev_of_unnormed_fitness_per_construct_related_df = _generate_stddev_of_unnormed_fitness_per_construct_related_df(
        fitted_per_rep_data_list)
    stddev_unnormed_fitness_per_construct_series = stddev_of_unnormed_fitness_per_construct_related_df["stddev_of_fc"]

    unnormed_fitness_posterior_prob_per_construct_series = _generate_posterior_prob_of_fitness_per_construct(
        unnormed_fitness_per_construct_series, stddev_of_unnormed_fitness_per_construct_related_df)

    # i.e., allbad, fc, sdfc, pp_fc
    return constructs_to_ignore_across_replicates_series, unnormed_fitness_per_construct_series, \
           stddev_unnormed_fitness_per_construct_series, unnormed_fitness_posterior_prob_per_construct_series


# Function generalizes functionality from Roman's analyze_dual-crispr_NA_combined_simple-null-w-lfdr.R line 35,
# shown below with additional comments by ab:
#    # ab: allbad is true for all constructs that lack at least 2 acceptable-abundance timepoints in BOTH experiments
#    allbad<-apply(good1,1,sum)<2 & apply(good2,1,sum)<2 #in this case I have nothing to use in either experiment
# Note that the R code requires 2 and only 2 replicates, while this code is agnostic as to replicate number.
def _get_constructs_to_ignore_across_replicates_mask(thresholded_per_rep_data_list):
    constructs_to_ignore_across_replicates_series = None
    for curr_per_replicate_info in thresholded_per_rep_data_list:
        if constructs_to_ignore_across_replicates_series is None:
            constructs_to_ignore_across_replicates_series = \
                curr_per_replicate_info._fails_min_num_timepts_bool_by_construct_series
        else:
            constructs_to_ignore_across_replicates_series = \
                constructs_to_ignore_across_replicates_series & \
                curr_per_replicate_info._fails_min_num_timepts_bool_by_construct_series

    return constructs_to_ignore_across_replicates_series  # i.e., allbad


# Function generalizes functionality from Roman's analyze_dual-crispr_NA_combined_simple-null-w-lfdr.R line 66,
# shown below with additional comments by ab:
#       # ab: i = construct index
#       # ab: v1 = variance of timepoints (e.g. mean number of days) for all the good timepoints for this construct in
#       # ab: replicate 1, v2 is analogous
#       # ab: f1 = covariance of log2 frequencies for all the good timepoints for this construct in replicate 1 with
#       # ab: the timepoints for those good timepoints, f2 is analogous
#       # ab: fc[i] is the combined fitness (across replicates) for construct i
#       fc[i]<-(f1+f2)/(v1+v2) #the combined fitness from replicate 1+2
#       #fc remains defined up to an additive constant
# Note that the R code requires 2 and only 2 replicates, while this code is agnostic as to replicate number, and that
# R code performs operations for a single construct at a time, while the code here does calculations for all constructs
# at the same time.
def _generate_unnormed_fitness_per_construct_series(thresholded_per_rep_data_list):
    descriptive_stats_per_construct_by_rep_num = _make_multiindex_df_across_replicates(
        thresholded_per_rep_data_list, "descriptive_stats_df")

    sum_covariances_per_construct_across_replicates_series = _calc_sum_of_stat_per_construct_across_replicates_series(
        descriptive_stats_per_construct_by_rep_num,
        ns_thresholded.ThresholdedPerReplicateData.get_covariance_of_log2_fractions_w_timepoints_header())
    sum_variances_per_construct_across_replicates_series = _calc_sum_of_stat_per_construct_across_replicates_series(
        descriptive_stats_per_construct_by_rep_num,
        ns_thresholded.ThresholdedPerReplicateData.get_variance_of_usable_timepoints_header())

    fitness_per_construct_across_replicates_series = sum_covariances_per_construct_across_replicates_series.div(
        sum_variances_per_construct_across_replicates_series)
    fitness_per_construct_across_replicates_series.fillna(0, inplace=True)  # NaNs will become 0 per Roman's code
    return fitness_per_construct_across_replicates_series  # i.e., fc


# Function generalizes functionality from Roman's analyze_dual-crispr_NA_combined_simple-null-w-lfdr.R line 103-112,
# shown below with additional comments by ab:
#       # ab: i = construct index
#       # ab: g1 = true/false values of whether construct i passes various abundance filters for all
#       # ab: timepoints for the first replicate; g2 is analogous
#       g1 <- good1[i,]
#       g2 <- good2[i,]
#       # ab: df = vector w one entry for each construct, with value from 0 to -2; 0 if both replicates of relevant
#       # ab: construct are good, -1 if just one is, -2 if neither are
#       # ab: I suspect that "df" is "degrees of freedom" for each construct
#       df[i]<-sum(g1)+sum(g2)-2
#       # ab: reminder: sqrtsum<-function(y) sqrt(sum(y^2))
#       # ab: outside of this function, this sdfc value is used only in the output of the construct file
#       # ab: I think that sdfc is "standard deviation of fc" for each construct.
#
#       # ab: So, Roman says:
#       # ab: std err of fc = sqrt of sum over t of (lower-case-epsilon for construct c, as a function of t)^2
#       # ab: divided by sqrt of (nc - 2)*sum over t of (t^2 - (mean of t)^2)
#       # ab: Note that lower-case-epsilon is Xc(t) - xc(t) = xX - xfitX
#       # ab: so xfitX - xX = - lower-case-epsilon ... but since it is being squared and then square-rooted, I
#       # ab: suppose the negation doesn't matter.
#       # ab: nc - 2 is the number of degrees of freedom
#       # ab: where nc = number of data points = 2*nt minus any number of points below the threshold
#       # ab: (note the description above seems to assume 2 replicates) ... seems to me this must mean
#       # ab: number of data points *for this construct c*, not total.
#       # ab: nt in above is number of timepoints, as here ...
#       # ab: Roman also says that tc [i.e., the t statistic for construct c] = fc /SE(fc)
#       # ab: and the internet tells me that SE(x) = SD(x)/sqrt(n) ... but maybe the sqrt(n) is just the most usual
#       # ab: sqrt of degrees of freedom, and could be something else in a more complex system.
#       # ab: So, the value being calculated directly below is SD(x), which is why it doesn't have the
#       # ab: nc - 2 term that is in the denominator of the SE(x) calculation (in the manuscript, Roman says that
#       # ab: fc's sd = sqrt(nc-2)*SE(fc), where nc-2 = degrees of freedom.  Farther below, after the calculation
#       # ab: of sdfc, we get the calculation of tc, which is fc/(sdfc/sqrt(df)), where the sdfc/sqrt(df) term
#       # ab: is the calculation of SE(fc).
#
#       # ab: result is vector with one std dev of fc for each construct
# #      sefc[i]<-sqrtsum( c(xfit1[i,good1[i,]],xfit2[i,good2[i,]]) - c(x1[i,good1[i,]],x2[i,good2[i,]]) ) /sqrtsum( c(time[good1[i,]],time[good2[i,]]) - mean(c(time[good1[i,]],time[good2[i,]])) ) /sqrt(df[i])
#       sdfc[i]<-sqrtsum( c(xfit1[i,g1],xfit2[i,g2]) - c(x1[i,g1],x2[i,g2]) ) /sqrtsum( c(time[g1],time[g2]) - mean(c(time[g1],time[g2])) )
#    }
# #what is median sd?
#    has_sd<-df>0
#    median_sd<-median(sdfc[has_sd])
#    sdfc[!has_sd]<-median_sd #just so it isn't 0
#
# Note that the R code requires 2 and only 2 replicates, while this code is agnostic as to replicate number.
def _generate_stddev_of_unnormed_fitness_per_construct_related_df(fitted_per_rep_data_list):
    log2_fract_per_construct_per_timept_per_replicate_df = _make_multiindex_df_across_replicates(
        fitted_per_rep_data_list, "_log2_fractions_by_constructs_by_timepoints_df")
    expected_log2_fract_per_construct_per_timept_per_replicate_df = _make_multiindex_df_across_replicates(
        fitted_per_rep_data_list, "_expected_log2_fraction_per_construct_per_timepoint_df")
    usable_constructs_by_timepoint_mask_per_replicate = _make_multiindex_df_across_replicates(
        fitted_per_rep_data_list, "_usable_constructs_by_timepoint_mask")
    concatenated_timepoints_series_per_replicate = _make_multiindex_df_across_replicates(fitted_per_rep_data_list,
                                                                                          "timepoints_series", axis=0)

    # stddev_of_fc = numpy array len of num of constructs filled with 0.1
    num_constructs = len(log2_fract_per_construct_per_timept_per_replicate_df.index.values)
    degrees_freedom = numpy.zeros(num_constructs)
    stddev_of_fc = numpy.empty(num_constructs)
    stddev_of_fc.fill(0.1)

    # for each construct
    idx = pandas.IndexSlice
    for curr_index, curr_construct in enumerate(log2_fract_per_construct_per_timept_per_replicate_df.index.values):
        curr_construct_slice = pandas.IndexSlice[curr_construct]
        usable_timepoints_per_replicate_mask = usable_constructs_by_timepoint_mask_per_replicate.loc[
            curr_construct_slice, idx[:, :]]
        if usable_timepoints_per_replicate_mask.any():  # if any of the timepoints are usable (TRUE)
            good_samples_slice = pandas.IndexSlice[:, usable_timepoints_per_replicate_mask]

            # TODO: Check w/Roman that the "2" being subtracted here in R code is for number of replicates
            degrees_freedom[curr_index] = usable_timepoints_per_replicate_mask.sum() - len(fitted_per_rep_data_list)

            curr_expected_log2_fracts_for_good_replicates = \
                expected_log2_fract_per_construct_per_timept_per_replicate_df.loc[
                    curr_construct_slice, good_samples_slice]
            curr_log2_fracts_for_good_replicates = log2_fract_per_construct_per_timept_per_replicate_df.loc[
                curr_construct_slice, good_samples_slice]
            numerator_dif = curr_expected_log2_fracts_for_good_replicates - curr_log2_fracts_for_good_replicates
            numerator = _calculate_sqrt_of_sum_of_squares(numerator_dif)

            curr_timepts_for_good_replicates = concatenated_timepoints_series_per_replicate[
                usable_timepoints_per_replicate_mask]
            denominator_diff = curr_timepts_for_good_replicates - curr_timepts_for_good_replicates.mean()
            denominator = _calculate_sqrt_of_sum_of_squares(denominator_diff)

            stddev_of_fc[curr_index] = numerator / denominator
        # end if
    # next construct

    has_stddev = degrees_freedom > 0
    median_sd = numpy.median(stddev_of_fc[has_stddev])
    stddev_of_fc[numpy.logical_not(has_stddev)] = median_sd  # Comment from Roman: just so it isn't 0
    result = pandas.DataFrame({"degrees_freedom":degrees_freedom,
                               "has_stddev":has_stddev,
                               "stddev_of_fc":stddev_of_fc},
                              index=log2_fract_per_construct_per_timept_per_replicate_df.index)
    return result


# Function implements analyze_dual-crispr_NA_combined_simple-null-w-lfdr.R line 19:
# sqrtsum<-function(y) sqrt(sum(y^2))
def _calculate_sqrt_of_sum_of_squares(one_d_array):
    return math.sqrt(numpy.sum(numpy.square(one_d_array)))


# Function generalizes functionality from Roman's analyze_dual-crispr_NA_combined_simple-null-w-lfdr.R line 115-122
# as well as 442, shown below with additional comments by ab:
#    # ab: for 1 to number of constructs
#    for (i in 1:nx) {
#       # ab: don't try to calculate t statistic and p value for any fc that doesn't have a stderr
#       if (!has_sd[i]) next
#       # ab: calc t statistic of fc
#       tstat[i]<-fc[i]/(sdfc[i]/sqrt(df[i]))
#       # ab: pt is R function, presumably to do t-test :)  Suppose pt result must be one-tailed, hence * 2
#       p_t[i]<-2*pt(-abs(tstat[i]),df=df[i]) #raw p-values from t-test
#    }
#
#    # ab: ok, lfdr stands for "local fdr", and the lfdr function comes from the qvalue bioconductor package
#    # ab: that is installed way above.  The first input is the vector of p-values, and the second input is the
#    # ab: estimated proportion of true null p-values, where pi0.method is "the method for automatically choosing tuning
#    # ab: parameter in the estimation of π0, the proportion of true null hypotheses."
#    # ab: nx = number of constructs; default value of lfdr is set to one for all of them
#    lfdr_fc<-rep(1,nx)
#    # ab: for constructs that have an sd, the lfdr of the fc could be calculated and is now set to its calculated value
#    # ab: instead of the default.
#    # ab: lfdr function comes from library(qvalue)
#    l<-lfdr(p_t[has_sd],pi0.method="bootstrap")
#    lfdr_fc[has_sd]<-l
#
#    # ab: I believe this is posterior probability of fc of each construct c (calculated across both replicates)
#    pp_fc<-1-lfdr_fc
# Note that the R code calculates tstat and p_t for a single construct at a time, while the code here does calculations
# for all constructs at the same time.
def _generate_posterior_prob_of_fitness_per_construct(fc_series, stddev_of_unnormed_fitness_per_construct_related_df):
    has_stddev_per_construct = stddev_of_unnormed_fitness_per_construct_related_df["has_stddev"].values
    stddev_of_fc_per_construct = stddev_of_unnormed_fitness_per_construct_related_df["stddev_of_fc"].values
    degs_free_per_construct = stddev_of_unnormed_fitness_per_construct_related_df["degrees_freedom"].values

    tstat_per_construct = fc_series.values/ (stddev_of_fc_per_construct / numpy.sqrt(degs_free_per_construct))
    raw_p_value_per_construct = 2 * scipy.stats.t.sf(numpy.absolute(tstat_per_construct), degs_free_per_construct)

    raw_p_value_per_construct[numpy.logical_not(has_stddev_per_construct)] = 1
    # TypeChecker inspection suppressed for statement including len(raw_p_value_per_construct) because PyCharm thinks
    # raw_p_value_per_construct is an integer when in fact it is a numpy array because scipy.stats.t.sf uses
    # "broadcasting"
    # noinspection PyTypeChecker
    lfdr_fc = numpy.ones(len(raw_p_value_per_construct))  # default value of lfdr is set to one for all constructs

    l = _calculate_lfdr_in_r(raw_p_value_per_construct[has_stddev_per_construct])
    lfdr_fc[has_stddev_per_construct] = l
    pp_fc = 1 - lfdr_fc

    return pandas.Series(pp_fc, index=fc_series.index)


def _calculate_lfdr_in_r(raw_pvals_per_construct_for_has_stddev):
    # TODO: All paths in this function need to be cleaned up so is not hard-coded like this but rather the way config is
    expected_output_fp = "/Users/Birmingham/Work/Repositories/mali-dual-crispr-pipeline/calc_lfdr_from_pval_out.txt"

    input_text = _process_pval_array_to_string(raw_pvals_per_construct_for_has_stddev)
    f = tempfile.NamedTemporaryFile()
    f.write(bytes(input_text, 'UTF-8'))

    try:
        process = subprocess.Popen(["bash", "/Users/Birmingham/Work/Repositories/mali-dual-crispr-pipeline/calc_lfdr_from_pval.bat",
                                    f.name, expected_output_fp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stderr:
            raise Exception(stderr)

        lfdr_per_construct_w_stddev_df = pandas.read_csv(expected_output_fp, sep="\t", header=0)
    finally:
        f.close()  # Once tempfile is closed, it gets removed

        # can't pull the same trick as above to delete R output file as R doesn't have tempfiles that I know of ...
        with contextlib.suppress(FileNotFoundError):
            os.remove(expected_output_fp)

    return numpy.asarray(lfdr_per_construct_w_stddev_df.iloc[:, 0].values)


def _process_pval_array_to_string(raw_pvals_per_construct_for_has_stddev):
    numpy.set_printoptions(threshold=numpy.inf, linewidth=numpy.inf)  # turn off summarization, line-wrapping
    output_text = numpy.array2string(raw_pvals_per_construct_for_has_stddev, separator='\n',
                                     suppress_small=True, precision=10)
    output_text = output_text.replace("[", "")
    output_text = output_text.replace("]", "")
    output_text = output_text.replace(" ", "")
    return output_text


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
