# third-party libraries
import numpy
import pandas
import scipy.stats

# ccbb libraries
import dual_crispr.construct_file_extracter as ns_extracter
import dual_crispr.probe_and_probe_pair_fitness_estimator as ns_probe_score

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


# FYI: Throughout, "tuple" is used as shorthand to refer to a tuple of (gene_name, probe_name) that uniquely
# identifies a gene/probe pair.  Thus, a data object that is "per_tuple" is per probe (given the assumption that
# each probe can only be for one gene), while a data object that is "by_tuple_by_tuple" is matrix-like and symmetric,
# with each axis indexed with entries for all gene/probe pairs in the experiment, and cell values being values of
# whatever the measured data is for the construct made up of the intersecting gene/probe-gene/probe pair.  Note that not
# all possible gene/probe-gene/probe pairs exist in actual constructs, so some cells will have NaN or default (usually
# zero) values.
#
# Alternately, when I designate a data object "per_probe" or "by_probe_by_probe", I mean that it is NOT indexed--these
# objects are generally either one-d or two-d numpy arrays, respectively. The dimensions of their axes represent the
# exact same thing as the "per_tuple" or "by_tuple_by_tuple" data objects (and in the same order), but since they are
# numpy arrays the slots along those dimensions are not explicitly labeled with what gene/probe pair they pertain to.
# Thus, these should be used with care.


def temp_do_irls_fitting(constructs_to_ignore_across_replicates_series, unnormed_fitness_per_construct_series,
                         stddev_unnormed_fitness_per_construct_series,
                         unnormed_fitness_posterior_prob_per_construct_series, scoring_data, num_iterations,
                         use_seed=False):
    # from temp_fit_ac_fc, pivoted and made symmetric
    unnormed_construct_fitnesses_by_tuple_by_tuple_df = scoring_data.make_by_tuple_by_tuple_annotated_symmetric_df(
        unnormed_fitness_per_construct_series)
    unnormed_fitness_posterior_prob_by_tuple_by_tuple_df = scoring_data.make_by_tuple_by_tuple_annotated_symmetric_df(
        unnormed_fitness_posterior_prob_per_construct_series)
    unnormed_fitness_sttdev_by_tuple_by_tuple_df = scoring_data.make_by_tuple_by_tuple_annotated_symmetric_df(
        stddev_unnormed_fitness_per_construct_series)
    ignore_mask_by_tuple_by_tuple_df = scoring_data.make_by_tuple_by_tuple_annotated_symmetric_df(
        constructs_to_ignore_across_replicates_series)
    is_nontargeting_mask_per_tuple_series = scoring_data.make_nontargeting_mask_per_tuple_series()

    initial_weights_by_tuple_by_tuple_df = _get_initial_weights_by_tuple_by_tuple_df(ignore_mask_by_tuple_by_tuple_df)

    mean_non_targeting_probe_fitness, normed_fitness_per_tuple_series, pi_scores_by_tuple_by_tuple_df = \
        _norm_fitness_and_calc_pi_scores(unnormed_construct_fitnesses_by_tuple_by_tuple_df,
                                         initial_weights_by_tuple_by_tuple_df,
                                         is_nontargeting_mask_per_tuple_series)

    # TODO: This would be incorrect for single-probe construct case where fitness_per_construct IS fitness_per_probe.
    # In that case, we just need normed_fitness_per_construct_series to be set equal to normed_fitness_per_tuple_series
    # (because normed_fitness_per_construct_series is returned but normed_fitness_per_tuple_series isn't).
    # In dual-crispr, times TWO because each construct has TWO probes in it, so subtract out effects of 2 null probes
    normed_fitness_per_construct_series = unnormed_fitness_per_construct_series - \
                                          2 * mean_non_targeting_probe_fitness

    rank_per_tuple_series = _rank_probes_within_genes(normed_fitness_per_tuple_series)
    mean_weighted_probe_fitness_per_gene_series = _calc_mean_weighted_probe_fitness_per_gene_series(
        normed_fitness_per_tuple_series, rank_per_tuple_series)

    # TODO: Find out from Roman *why* mean_weighted_pi_stacked_by_probe_pair_df is calculated here when never used?
    # TODO: assume these calls would be left out in single-probe construct case.
    rank_weightings_by_tuple_by_tuple_for_expressed_df = _get_rank_weightings_by_tuple_by_tuple_for_expressed_df(
        rank_per_tuple_series, initial_weights_by_tuple_by_tuple_df)
    mean_weighted_pi_stacked_by_probe_pair_df = _calc_mean_weighted_pi_stacked_by_probe_pair_df(
        pi_scores_by_tuple_by_tuple_df, rank_weightings_by_tuple_by_tuple_for_expressed_df)

    mean_weighted_pi_results_list = []
    mean_weighted_probe_fitness_per_gene_results_list = []
    for curr_iteration in range(num_iterations):
        # TODO: recreate status print-out from Roman's code?
        # cat("\n", iter, "\n")

        perturbed_construct_fitnesses_by_tuple_by_tuple_df = _perturb_fitness_by_tuple_by_tuple_df(
            unnormed_construct_fitnesses_by_tuple_by_tuple_df, unnormed_fitness_posterior_prob_by_tuple_by_tuple_df,
            unnormed_fitness_sttdev_by_tuple_by_tuple_df)

        # TODO: Again, think that in single-probe construct case, the only thing we need to do here is normalize
        # perturbed_construct_fitnesses_by_tuple_by_tuple_df using non-targeting probe info, because per-construct
        # fitness IS per-probe fitness
        curr_mean_non_targeting_probe_fitness, curr_normed_fitness_per_tuple_series, \
        curr_pi_scores_by_tuple_by_tuple_df = _norm_fitness_and_calc_pi_scores(
            perturbed_construct_fitnesses_by_tuple_by_tuple_df, initial_weights_by_tuple_by_tuple_df,
            is_nontargeting_mask_per_tuple_series)

        curr_mean_weighted_probe_fitness_per_gene_series = _calc_mean_weighted_probe_fitness_per_gene_series(
            curr_normed_fitness_per_tuple_series, rank_per_tuple_series)
        mean_weighted_probe_fitness_per_gene_results_list.append(curr_mean_weighted_probe_fitness_per_gene_series)

        # TODO: I think that these calls would be left out in single-probe construct case; check with Roman
        curr_mean_weighted_pi_stacked_by_probe_pair_df = _calc_mean_weighted_pi_stacked_by_probe_pair_df(
            curr_pi_scores_by_tuple_by_tuple_df, rank_weightings_by_tuple_by_tuple_for_expressed_df)
        mean_weighted_pi_results_list.append(curr_mean_weighted_pi_stacked_by_probe_pair_df)

    pi_null, f_sd, gene_pair_statistics_df = _calc_stats_from_combined_iteration_info(
        mean_weighted_probe_fitness_per_gene_series, mean_weighted_probe_fitness_per_gene_results_list,
        mean_weighted_pi_results_list)

    # TODO: In single-probe construct case, pi_null would be irrelevant; how to handle?
    # also, in that case normed_fitness_per_construct_series is in fact equal to normed_fitness_per_tuple_series.
    return pi_null, mean_weighted_probe_fitness_per_gene_series, rank_per_tuple_series, f_sd, \
           normed_fitness_per_construct_series, gene_pair_statistics_df


def _get_initial_weights_by_tuple_by_tuple_df(ignore_mask_by_tuple_by_tuple_df):
    ignore_mask_by_probe_by_probe_array = ignore_mask_by_tuple_by_tuple_df.values
    # 0 for all bad/non-existent probe intersections (constructs) and 1 for all usable constructs
    initial_weights_by_probe_by_probe_array = numpy.logical_not(ignore_mask_by_probe_by_probe_array).astype(numpy.int)
    initial_weights_by_tuple_by_tuple_df = pandas.DataFrame(initial_weights_by_probe_by_probe_array,
                                                            index=ignore_mask_by_tuple_by_tuple_df.index,
                                                            columns=ignore_mask_by_tuple_by_tuple_df.columns)
    return initial_weights_by_tuple_by_tuple_df


# a
# TODO: I think that in single-probe construct case, the call to _calc_unnormed_probe_fitnesses_and_probe_pair_pi_scores
# would be irrelevant: the unnormed_construct_fitnesses_by_tuple_by_tuple_df that come in (assuming they are
# really unnormed_construct_fitnesses_by_tuple_df, a series represented as a df) ARE in fact the
# unnormed_fitness_per_tuple_series, so only the norming step would be necessary.
def _norm_fitness_and_calc_pi_scores(unnormed_construct_fitnesses_by_tuple_by_tuple_df,
                                     initial_weights_by_tuple_by_tuple_df,
                                     is_nontargeting_mask_per_tuple_series):

    unnormed_fitness_per_tuple_series, pi_scores_by_tuple_by_tuple_df = \
        ns_probe_score._calc_unnormed_probe_fitnesses_and_probe_pair_pi_scores(
        unnormed_construct_fitnesses_by_tuple_by_tuple_df, initial_weights_by_tuple_by_tuple_df)

    # Get unnormed_fitness for non-targeting probes only and take mean of their fitnesses (this is a scalar)
    mean_non_targeting_probe_fitness = unnormed_fitness_per_tuple_series[is_nontargeting_mask_per_tuple_series].mean()
    # real, normalized fitness of any probe is the fitness we see minus the fitness we expect for probes that do nothing
    # (i.e., non-targeting probes)
    normed_fitness_per_tuple_series = unnormed_fitness_per_tuple_series - mean_non_targeting_probe_fitness

    return mean_non_targeting_probe_fitness, normed_fitness_per_tuple_series, pi_scores_by_tuple_by_tuple_df


# b
def _rank_probes_within_genes(normed_fitness_per_tuple_series):
    normed_fitness_per_gene_and_probe_series = _transform_per_tuple_to_gene_and_probe_hier_index_series(
        normed_fitness_per_tuple_series)
    # TODO: refactor way of specifying level name
    normed_fitness_per_probe_grouped_by_gene = normed_fitness_per_gene_and_probe_series.groupby(level=['gene'])
    rank_per_gene_and_probe_series = normed_fitness_per_probe_grouped_by_gene.transform(_rank_probes_for_a_gene)
    rank_per_tuple_series = _transform_per_gene_and_probe_hierarchical_to_per_tuple_indexing_series(
        rank_per_gene_and_probe_series)
    return rank_per_tuple_series


def _rank_probes_for_a_gene(a_series):
    # TODO: refactor code to check if the gene for which we're getting data is actually the "nontargeting" gene--
    # probably need to ask scoring_data
    multiplier = -1 if a_series.name == "nontargeting" else 1

    rankable_series = multiplier * a_series.abs()
    # PyCharm thinks rankable_series is an int, not a series, thus the need to suppress inspection that gives
    # incorrect error
    # noinspection PyUnresolvedReferences
    ascending_ranks = scipy.stats.rankdata(rankable_series.values, method="ordinal")
    descending_zero_based_ranks = len(rankable_series.values) - ascending_ranks
    return pandas.Series(descending_zero_based_ranks, index=a_series.index)


# c
def _get_rank_weightings_by_tuple_by_tuple_for_expressed_df(rank_per_tuple_series,
                                                            initial_weights_by_tuple_by_tuple_df):
    rank_per_tuple_df = rank_per_tuple_series.to_frame()
    rank_weightings_by_tuple_by_tuple_df = rank_per_tuple_df.dot(rank_per_tuple_df)
    rank_weightings_by_tuple_by_tuple_for_expressed_df = rank_weightings_by_tuple_by_tuple_df[
        initial_weights_by_tuple_by_tuple_df > 0]
    return rank_weightings_by_tuple_by_tuple_for_expressed_df


# d
# think this method should work fine with single-probe constructs
def _calc_mean_weighted_probe_fitness_per_gene_series(normed_fitness_per_tuple_series,
                                                      rank_per_tuple_series):

    # TODO: refactor way of specifying new column names throughout
    # TODO: can I actually use numpy.square on a series (and get BACK a series?)
    weighted_fitness_per_probe_series = numpy.square(rank_per_tuple_series) * normed_fitness_per_tuple_series
    normed_fitness_and_rank_per_tuple_df = pandas.concat([weighted_fitness_per_probe_series, rank_per_tuple_series],
                                                             axis=1)
    normed_fitness_and_rank_per_tuple_df.columns = ["fitness", "rank"]

    normed_fitness_and_rank_per_gene_and_probe_df = _transform_per_tuple_to_gene_and_probe_hier_index_df(
        normed_fitness_and_rank_per_tuple_df)
    # TODO: refactor way of specifying level name
    normed_fitness_and_rank_per_probe_grouped_by_gene = normed_fitness_and_rank_per_gene_and_probe_df.groupby(
        level="gene")
    mean_weighted_probe_fitness_per_gene_series = (normed_fitness_and_rank_per_probe_grouped_by_gene["fitness"].sum() /
                                                   normed_fitness_and_rank_per_probe_grouped_by_gene["rank"].sum())
    return mean_weighted_probe_fitness_per_gene_series


# e
# TODO: this method needs more work to nail down the logic
def _calc_mean_weighted_pi_stacked_by_probe_pair_df(pi_scores_by_tuple_by_tuple_df,
                                                    rank_weightings_by_tuple_by_tuple_for_expressed_df,
                                                    small_positive_value=1e-6):
    raise NotImplementedError

    # multiply the two dataframes element-wise
    weighted_pi_by_tuple_by_tuple_for_expressed_df = pi_scores_by_tuple_by_tuple_df * \
                                                     rank_weightings_by_tuple_by_tuple_for_expressed_df

    # TODO: really need to group by gene in *both* dimensions and to sum *everything* in the little square
    # TODO: re-visit this logic, as I am not sure it is right.  The code below is intended to fill the by_gene_by_gene
    # matrix for *every* pair of genes, but Roman's code sets default by_gene_by_gene matrix values to zero and then
    # fills in values for ONLY the upper triangle (with diagonal left as zero).  Kind of confusing, because later the
    # code *takes* the upper triangle of this matrix that only has the upper triangle filled in anyway ...
    # I get as a result.  How to make that happen?
    numerator = weighted_pi_by_probe_by_probe_for_expressed_df.groupby(level="gene").sum()
    denominator = rank_weightings_by_tuple_by_tuple_for_expressed_df.groupby(level="gene").sum()
    denominator[denominator < small_positive_value] = small_positive_value

    mean_weighted_pi_by_gene_by_gene_df = numerator / denominator
    # drop row and column for non-targeting "gene"
    # TODO refactor way of specifying nontargeting gene
    mean_weighted_pi_by_gene_by_gene_df.drop("nontargeting", axis=0)  # drop row
    mean_weighted_pi_by_gene_by_gene_df.drop("nontargeting", axis=1)  # drop column

    by_gene_by_gene_upper_tri_mask = _get_upper_tri_wo_diagonal(numpy.ones(
        mean_weighted_pi_by_gene_by_gene_df.shape).astype(numpy.bool))
    mean_weighted_pi_by_gene_by_gene_targeting_only_upper_tri_df = mean_weighted_pi_by_gene_by_gene_df.where(
        by_gene_by_gene_upper_tri_mask)
    mean_weighted_pi_stacked_by_probe_pair_df = mean_weighted_pi_by_gene_by_gene_targeting_only_upper_tri_df.stack(). \
        reset_index()
    # TODO: refactor assigned column names
    mean_weighted_pi_stacked_by_probe_pair_df.columns = ['probeA', 'probeB', 'Value']
    mean_weighted_pi_stacked_by_probe_pair_df.set_index(['probeA', 'probeB'], inplace=True)
    return mean_weighted_pi_stacked_by_probe_pair_df


def _get_upper_tri_wo_diagonal(matrix_df):
    # k = 1 means main diagonal is NOT included
    return pandas.DataFrame(numpy.triu(matrix_df, k=1), index=matrix_df.index, columns=matrix_df.columns)


def _make_symmetric_matrix_df_from_upper_tri_wo_diag_df(upper_tri_wo_diag_df, diagonals_value=0):
    # note that fill_value is only used if at least one of the values is NOT missing; if both are missing, the
    # result is na
    result_df = upper_tri_wo_diag_df.add(upper_tri_wo_diag_df.transpose(), fill_value=0)
    numpy.fill_diagonal(result_df.values, diagonals_value)
    return result_df


def _perturb_fitness_by_tuple_by_tuple_df(unnormed_construct_fitnesses_by_tuple_by_tuple_df,
                                          unnormed_fitness_posterior_prob_by_tuple_by_tuple_df,
                                          stddev_of_fitness_by_probe_by_probe_df):
    # TODO: write code to perturb construct_fitnesses_by_tuple_by_tuple_df:
    # same starting value as fc_0<-matrix(0,nrow=nprobes,ncol=nprobes)
    # fc_1 <- matrix(0, nrow = nprobes, ncol = nprobes)
    #
    # # set seed ONLY FOR TESTING
    # if (useSeed == TRUE) {
    #   set.seed(iter)
    # }
    # # ok, previous fc0 was fc0<-fc ; now we're adding a random normal
    # # to each fc, where each normal variable's mean is ?the sum of all
    # # fcs in the upper triangle?, and its stddev is the stddev of the
    # # fc of the analogous construct?
    # fc0 <- fc_0[utri] + rnorm(ntri, sd = sdfc_0[utri])
    # # gets the posterior probability of the fcs for constructs in the
    # # upper triangle
    # pp0 <- pp_0[utri]
    #
    # # set seed ONLY FOR TESTING
    # if (useSeed == TRUE) {
    #   set.seed(iter)
    # }
    #
    # # draw is 1 if the random posterior probability is less than the
    # # calculated posterior probability, zero otherwise
    # draw <- ifelse(runif(ntri) < pp0, 1, 0)
    # # ok, multiplying by draw will set to zero every value in fc_1
    # # where the posterior probability of the real fc was not more
    # # than you'd expect by chance
    # fc_1[utri] <- fc0 * draw
    # fc_1 <- fc_1 + t(fc_1) # make fc_1 symmetric

    raise NotImplementedError
    return perturbed_construct_fitnesses_by_tuple_by_tuple_df


def _transform_per_tuple_to_gene_and_probe_hier_index_series(per_tuple_series):
    # TODO: write code to change from per_tuple indexing to hierarchical gene/probe indexing so can group by gene
    raise NotImplementedError
    return per_gene_and_probe_series


def _transform_per_tuple_to_gene_and_probe_hier_index_df(per_tuple_df):
    # TODO: write code to change from per_tuple indexing to hierarchical gene/probe indexing so can group by gene
    raise NotImplementedError
    return per_gene_and_probe_df


def _transform_per_gene_and_probe_hierarchical_to_per_tuple_indexing_series(per_gene_and_probe_series):
    # TODO: write code to change back to per_tuple indexing
    raise NotImplementedError
    return per_tuple_series


def _transform_by_tuple_by_tuple_to_by_gene_and_probe_by_gene_and_probe_hier_index(by_tuple_by_tuple_df):
    # TODO: write code to change from by_tuple_by_tuple indexing to hierarchical gene/probe by gene/probe indexing
    # so can group by gene
    raise NotImplementedError
    return by_gene_and_probe_by_gene_and_probe_df


# TODO: I think the only part of this function that applies in the single-gene construct case is the calculation of f_sd
def _calc_stats_from_combined_iteration_info(mean_weighted_probe_fitness_per_gene_series,
                                             mean_weighted_probe_fitness_per_gene_results_list,
                                             mean_weighted_pi_results_list):
    raise NotImplementedError
    # f = mean_weighted_probe_fitness_per_gene_series--e.g, single gene fitness series.  Note that this series comes
    # from the calculation done PREVIOUS TO and OUTSIDE the iteration loop: the
    # curr_mean_weighted_probe_fitness_per_gene_series created within the iteration loop are appended to
    # mean_weighted_probe_fitness_per_gene_results_list and used ONLY for the calculation of f_sd

    # f_iter = mean_weighted_probe_fitness_per_gene_results_list
    # pi_iter = mean_weighted_pi_results_list

    # "1" in apply call means apply the function (e.g., sd) over each row
    # f_sd <- apply(f_iter, 1, sd)  # output in single gene fitness file
    #
    # used in several plots, output in pi score file
    # pi_mean <- apply(pi_iter, 1, mean)
    #
    # pi_sd < - apply(pi_iter, 1, sd)  # output in pi score file
    #
    # used directly below, then *redefined* (same way) and used (differently) in prepping for pi score output file
    # pi_iter_null <- pi_iter - pi_mean
    #
    # used directly below, and in histogram of pi scores
    # pi_null <- c(pi_iter_null, -pi_iter_null)
    #
    # enull < - ecdf(pi_null)
    # emean < - ecdf(pi_mean)
    #
    # fdr_left < - pmin(1, enull(pi_mean) / emean(pi_mean))
    # fdr_right < - pmin(1, (enull(-pi_mean)) / (1 - emean(pi_mean)))
    #

    # genePairInfoDf = gatherGenePairInfo(preppedAndCheckedData$genes, pi_iter, uutri, pi_mean, pi_sd, fdr_left, fdr_right, f)
    # gatherGenePairInfo will become _generate_gene_pair_statistics_df . Note that there will be other input parameters
    # to _generate_gene_pair_statistics_df as well, but all the other non-trivial ones are generated within this calling
    #  function
    gene_pair_statistics_df = _generate_gene_pair_statistics_df(mean_weighted_probe_fitness_per_gene_series)
    return pi_null, f_sd, gene_pair_statistics_df


def _calc_ecdf(a_value):
    # TODO: Determine whether this ecdf calc to be done in R natively or can be done with scipy.stats or in R via rpy2
    raise NotImplementedError
    return ecdf_of_value


# I don't think any of this function applies in the single-probe construct case as there ARE no gene pairs.
def _generate_gene_pair_statistics_df(mean_weighted_probe_fitness_per_gene_series):
    # there will be other input parameters as well, but all the other non-trivial ones are generated in
    # _calc_stats_from_combined_iteration_info

    # f = mean_weighted_probe_fitness_per_gene_series--e.g, single gene fitness series
    # gatherGenePairInfo < - function(genes, pi_iter, uutri, pi_mean, pi_sd, fdr_left, fdr_right, f)
    # {
    #     n = length(genes)
    #
    # # So ... we're building 3 matrices of n genes by n genes.
    # # The first one just
    # # Again, here 1 & 2 hardcoding acceptable as means 1st and 2nd of DUAL crispr construct
    # g1names < - matrix("", ncol=n, nrow=n)  # n = number of genes
    # g2names < - matrix("", ncol=n, nrow=n)
    # ggnames < - matrix("", ncol=n, nrow=n)
    # for (i in 1: (n - 1)) {
    # for (j in (i + 1): n) {
    #     g1names[i, j] < - genes[i]
    # g2names[i, j] < - genes[j]
    # # TODO: Separator should be set in params rather than hardcoded
    # ggnames[i, j] < - paste(genes[i], "_", genes[j], sep="")
    # ggnames[j, i] < - ggnames[i, j]
    # }
    # }
    # # end make gene name matrices
    #
    # z < - pi_mean / sd(pi_mean)
    #
    # pi_iter_null < - pi_iter - pi_mean
    # abspi < - abs(pi_mean)
    # PP < - apply(abs(pi_iter_null) < abspi, 1, mean)
    #
    # # get names of these gene pairs
    # # These g1//g2, etc hardcodes are ok because represent DUAL crispr genes
    # names_of_g1 < - g1names[uutri]
    # fg1 < - f[names_of_g1]
    # names_of_g2 < - g2names[uutri]
    # fg2 < - f[names_of_g2]
    # fg12 < - fg1 + fg2
    # names_of_gg < - ggnames[uutri]
    #
    # res < - data.frame(
    #     names_of_gg=names_of_gg,
    #     names_of_g1=names_of_g1,
    #     fg1=fg1,
    #     names_of_g2=names_of_g2,
    #     fg2=fg2,
    #     fg12=fg12,
    #     pi_mean=pi_mean,
    #     pi_sd=pi_sd,
    #     PP=PP,
    #     abspi=abspi,
    #     fdr_left=fdr_left,
    #     fdr_right=fdr_right,
    #     z=z,
    #     row.names = names_of_gg
    # )
    #
    # return (res)
    raise NotImplementedError
    return gene_pair_statistics_df
