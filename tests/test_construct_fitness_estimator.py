# standard libraries
import unittest

# third-party libraries
import numpy
import pandas
import pandas.util.testing
import scipy.stats

# project-specific libraries
import dual_crispr.construct_fitness_estimator as ns_test
import dual_crispr.gene_and_gene_pair_fitness_estimator as ns_test_2
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

    def test_generate_unnormed_fitness_per_construct_series(self):
        input_thresh_rep_list = [ns_thresh.TestThresholdedPerReplicateData.help_make_thresh_rep(True),
                                 ns_thresh.TestThresholdedPerReplicateData.help_make_thresh_rep(False)]
        expected_output_series = ns_help.TestScorer.help_make_naive_fcs_series()

        real_output_series = ns_test._generate_unnormed_fitness_per_construct_series(input_thresh_rep_list)
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
        expected_output_df = ns_help.TestScorer.help_make_multiindex_df()
        real_output_df = ns_test._make_multiindex_df_across_replicates(input_thresh_rep_list, "descriptive_stats_df")

        # belt and suspenders approach: check values within a precision with code below, then check indexes, etc,
        # with assert_frame_equal below (with unrealistically strong rounding to make values match)
        ns_help.help_test_df_equality(self, expected_output_df, real_output_df, 0.0000005)

    def test__calc_sum_of_stat_per_construct_across_replicates_series(self):
        input_df = ns_help.TestScorer.help_make_multiindex_df()
        expected_output_series = ns_help.access_test_file_as_series("sum_covariance_of_log2_fractions_w_timepoints.txt",
                                                                    set_index=True)
        real_output_series = ns_test._calc_sum_of_stat_per_construct_across_replicates_series(
            input_df, "covariance_of_log2_fractions_w_timepoints")

        # belt and suspenders approach: check values within a precision with code below, then check indexes, etc,
        # with assert_frame_equal below (with unrealistically strong rounding to make values match)
        ns_help.help_test_series_equality(self, expected_output_series, real_output_series, 0.0000000005)

    def test__generate_stddev_of_unnormed_fitness_per_construct_related_df(self):
        fitted_per_rep_list = [ns_fitted.TestFittedPerReplicateData.help_make_fitted_rep(get_rep_1=True),
                               ns_fitted.TestFittedPerReplicateData.help_make_fitted_rep(get_rep_1=False)]
        expected_output_df = ns_help.TestScorer.help_make_stddev_of_unnormed_fitness_per_construct_related_df()

        real_output_df = ns_test._generate_stddev_of_unnormed_fitness_per_construct_related_df(fitted_per_rep_list)
        ns_help.help_test_df_equality(self, expected_output_df, real_output_df, 0.000001)


    def test__calculate_lfdr_in_r(self):
        expected_input_series = ns_help.access_test_file_as_series(
            "/Users/Birmingham/Work/Repositories/mali-dual-crispr-pipeline/dual_crispr/distributed_files/test_data/"
            "test_set_8/temp_p_t_for_has_stddev_from_R.txt", col_index_to_use_as_index=None, set_index=False)
        expected_output_series = ns_help.access_test_file_as_series("test_l.txt", col_index_to_use_as_index=None,
                                                                    set_index=False)
        real_output = ns_test._calculate_lfdr_in_r(expected_input_series.values)
        real_output_series = pandas.Series(real_output)
        pandas.util.testing.assert_series_equal(expected_output_series, real_output_series)

    def test__generate_fitness_raw_p_values(self):
        input_naive_fc_series = ns_help.TestScorer.help_make_naive_fcs_series()
        input_stddev_related_df = ns_help.TestScorer.help_make_stddev_of_unnormed_fitness_per_construct_related_df()
        expected_output_series = ns_help.access_test_file_as_series("6_postposteriorprobplot_pp_fc_new.txt",
                                                                    clear_index_name=True)

        real_output_series = ns_test._generate_posterior_prob_of_fitness_per_construct(input_naive_fc_series,
                                                                                       input_stddev_related_df)
        pandas.util.testing.assert_series_equal(expected_output_series, real_output_series)

    def test_pivot(self):
        construct_ones = ["probeA", "probeB", "probeC"]
        construct_twos = ["probeB", "probeC", "probeA"]
        fc_vals = [12.5, 18.39, 9.83]
        test_df = pandas.DataFrame({"probe_one": construct_ones, "probe_two": construct_twos, "fc": fc_vals})
        # pivot_table = test_df.pivot(index='probe_one', columns='probe_two', values='fc')
        pivot_table = ns_test_2._make_by_probe_by_probe_symmetric_df(test_df)
        self.fail("test is temporary")

    def test_group(self):
        arrays = [['brca1', 'brca1', 'brca1', 'nontargeting', 'nontargeting', 'nontargeting', 'sed', 'sed'],
                  ['brca1_1', 'brca1_2', 'brca1_3', 'nontargeting_1', 'nontargeting_2', 'nontargeting_3', 'sed_1', 'sed_2']]
        tuples = list(zip(*arrays))
        index = pandas.MultiIndex.from_tuples(tuples, names=['gene', 'probe'])
        s = pandas.Series([12,-15,19,10,135,-34,1643,45], index=index)

        expected_ranks_series = pandas.Series([2,1,0,0,2,1,0,1], index=index)

        def rank_values_ordinally(a_series):
            multiplier = -1 if a_series.name == "nontargeting" else 1
            rankable_series = multiplier * numpy.absolute(a_series)
            ascending_ranks = scipy.stats.rankdata(rankable_series.values, method="ordinal")
            descending_zero_based_ranks = len(rankable_series.values) - ascending_ranks
            return pandas.Series(descending_zero_based_ranks, index=a_series.index)

        s_by_gene = s.groupby(level=['gene'])
        rank_by_gene = s_by_gene.transform(rank_values_ordinally)
        pandas.util.testing.assert_series_equal(expected_ranks_series, rank_by_gene)
        #self.fail("test is temporary")

    def test_indexes(self):
        test_square_df = pandas.DataFrame({'a': [12, -15, 19, 10],
                                           'b': [1, 135, -34, 1643],
                                           'c': [1, 1, 45, 309],
                                           'd': [1, 1, 1, 92]}, index=['a', 'b', 'c', 'd'])

        test_square_2_df = pandas.DataFrame({'a': [2, -5, 9, 3],
                                           'b': [1, 35, -4, 643],
                                           'c': [1, 1, 5, 9],
                                           'd': [1, 1, 1, 2]}, index=['a', 'b', 'c', 'd'])

        upper_tri_mask = numpy.triu(numpy.ones(test_square_df.shape)).astype(numpy.bool)

        def reformat(test_square_df, upper_tri_mask):
            mod_test_df = test_square_df.where(upper_tri_mask)
            stacked_test_df = mod_test_df.stack().reset_index()
            stacked_test_df.columns = ['probeA', 'probeB', 'Value']
            indexed_test_df = stacked_test_df.set_index(['probeA', 'probeB'])
            return indexed_test_df

        reformatted_test_df = reformat(test_square_df, upper_tri_mask)
        reformatted_test_2_df = reformat(test_square_2_df, upper_tri_mask)

        reformatted_test_df["next_val"] = reformatted_test_2_df["Value"]
        self.fail("test is temporary")

    def test_multid_groupby(self):
        arrays = [['brca1', 'brca1', 'sed', 'sed', 'sept1', 'sept1'],
                  ['a', 'b', 'c', 'd', 'e', 'f']]
        tuples = list(zip(*arrays))
        index = pandas.MultiIndex.from_tuples(tuples, names=['gene', 'probe'])

        test_square_df = pandas.DataFrame([[12, -15, 19, 10, 82, 2],
                                            [1, 135, -34, 1643, 4, 12],
                                            [1, 1, 45, 309, 489, 23],
                                            [1, 1, 1, 92, 438, 2],
                                            [1, 1, 1, 1, 7, 23],
                                            [1, 1, 1, 1, 1, 9]], index=index, columns=index)
        grouped = test_square_df.groupby(level="gene", axis=0)

        for name, group in grouped:
            relevant_df = group.loc[pandas.IndexSlice[:, :], pandas.IndexSlice[name, :]]
            sum_all = relevant_df.values.sum()
            print(sum_all)
        self.fail("test is temporary")

    def test_index_setting(self):
        arrays = [['brca1', 'brca1', 'sed', 'sed', 'sept1', 'sept1'],
                  ['a', 'b', 'c', 'd', 'e', 'f']]
        tuples = list(zip(*arrays))
        index = pandas.MultiIndex.from_tuples(tuples, names=['gene', 'probe'])

        data_lists = [["brca1_a__sed_c", "brca1", "a", "sed", "c", -15],
                      ["brca1_b__sed_c", "brca1", "b", "sed", "c", 12],
                      ["brca1_a__sed_d", "brca1", "a", "sed", "d", 8],
                      ["brca1_b__sed_d", "brca1", "b", "sed", "d", 3],
                      ["sed_c__sept1_e", "sed", "c", "sept1", "e", 2],
                      ["sed_d__sept1_e", "sed", "d", "sept1", "e", 6],
                      ["sed_c__sept1_f", "sed", "c", "sept1", "f", 122],
                      ["sed_d__sept1_f", "sed", "d", "sept1", "f", 92],
                      ["brca1_a__sept1_e", "brca1", "a", "sept1", "e", -9],
                      ["brca1_b__sept1_e", "brca1", "b", "sept1", "e", 90],
                      ["brca1_a__sept1_f", "brca1", "a", "sept1", "f", 37],
                      ["brca1_b__sept1_f", "brca1", "b", "sept1", "f", 123]]

        data_lists2 = [["brca1_a__sed_c", -150],
                      ["brca1_b__sed_c", 120],
                      ["brca1_a__sed_d", 80],
                      ["brca1_b__sed_d", 30],
                      ["sed_c__sept1_e", 20],
                      ["sed_d__sept1_e", 60],
                      ["sed_c__sept1_f", 1220],
                      ["sed_d__sept1_f", 920],
                      ["brca1_a__sept1_e", -90],
                      ["brca1_b__sept1_e", 900],
                      ["brca1_a__sept1_f", 370],
                      ["brca1_b__sept1_f", 1230]]

        annotated_df = pandas.DataFrame(data_lists,
                                          columns=["construct", "geneA", "probeA", "geneB", "probeB", "values"])
        annotated_df.set_index(["construct"], inplace=True)
        to_annotate_df = pandas.DataFrame(data_lists2, columns=["construct", "fc"])
        to_annotate_df.set_index(["construct"], inplace=True)
        newly_annotated_df = pandas.concat([to_annotate_df, annotated_df[["geneA", "probeA", "geneB", "probeB"]]],
                                           axis=1)
        newly_annotated_df["a"] = list(zip(newly_annotated_df["geneA"], newly_annotated_df["probeA"]))
        newly_annotated_df["b"] = list(zip(newly_annotated_df["geneB"], newly_annotated_df["probeB"]))
        newly_annotated_pruned_df = newly_annotated_df.drop(["geneA", "probeA", "geneB", "probeB"], axis=1)
        pivoted_df = newly_annotated_pruned_df.pivot(index="a", columns="b", values="fc")
        symmetric_pivoted_df = pivoted_df.add(pivoted_df.transpose(), fill_value=0)
        symmetric_pivoted_df.fillna(0, inplace=True)
        numpy.fill_diagonal(symmetric_pivoted_df.values, 0)

        arrays = [newly_annotated_df.loc[:, "geneA"].values,
                  newly_annotated_df.loc[:, "probeA"].values]
        tuples = list(zip(*arrays))
        index = pandas.MultiIndex.from_tuples(tuples, names=["geneA", "probeA"])

        arrays2 = [newly_annotated_df.loc[:, "geneB"].values,
                  newly_annotated_df.loc[:, "probeB"].values]
        tuples2 = list(zip(*arrays2))
        index2 = pandas.MultiIndex.from_tuples(tuples2, names=["geneB", "probeB"])

        mod_annotated_df = newly_annotated_df[["geneB", "probeB", "fc"]]
        mod_annotated_df.index = index
        mod_annotated_df.columns = index2


        # grouped = test_square_df.groupby(level="gene", axis=0)
        #
        # for name, group in grouped:
        #     relevant_df = group.loc[pandas.IndexSlice[:, :], pandas.IndexSlice[name, :]]
        #     sum_all = relevant_df.values.sum()
        #     print(sum_all)
        self.fail("test is temporary")
