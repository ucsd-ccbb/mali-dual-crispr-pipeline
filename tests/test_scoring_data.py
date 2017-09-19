# standard libraries
import io
import unittest

# third-party libraries
import pandas
import pandas.util.testing

# project-specific libraries
import dual_crispr.per_replicate_data as ns_per_rep
import dual_crispr.scoring_data as ns_test
import tests.test_per_replicate_data as ns_help

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestScoringData(unittest.TestCase):
    @staticmethod
    def help_make_scoring_data(input_str):
        scoring_prepped_df = pandas.read_csv(io.StringIO(input_str), sep="\t")
        return ns_test.ScoringData(scoring_prepped_df, "NonTargeting", _run_init_methods=False)

    def test__rename_non_targeting_gene(self):
        input_str = """construct_id	target_a_id	probe_a_id	target_b_id	probe_b_id	target_pair_id	probe_pair_id
0	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	BRCA1_NonTargetingControlGuideForHuman0412	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412
1	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	NonTargetingControlGuideForHuman0352	SETD2	SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352_SETD2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972
2	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	BRCA1	BRCA1_chr17_41256141	NonTargetingControlGuideForHuman0362	NonTargetingControlGuideForHuman0362	BRCA1_NonTargetingControlGuideForHuman0362	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362
"""
        expected_output_str = """construct_id	target_a_id	probe_a_id	target_b_id	probe_b_id	target_pair_id	probe_pair_id
0	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	0	NonTargetingControlGuideForHuman0412	BRCA1_NonTargetingControlGuideForHuman0412	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412
1	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	0	NonTargetingControlGuideForHuman0352	SETD2	SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352_SETD2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972
2	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	BRCA1	BRCA1_chr17_41256141	0	NonTargetingControlGuideForHuman0362	BRCA1_NonTargetingControlGuideForHuman0362	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362
"""
        scoring_data = TestScoringData.help_make_scoring_data(input_str)
        scoring_data._rename_non_targeting_gene()

        expected_output = ns_help.help_make_df(expected_output_str)
        self.assertTrue(expected_output.equals(scoring_data._input_counts_and_annotation_df))

    def test__remove_constructs_w_both_probes_targeting_same_target_id(self):
        input_str = """construct_id	target_a_id	probe_a_id	target_b_id	probe_b_id	target_pair_id	probe_pair_id	A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
0	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	0	NonTargetingControlGuideForHuman0412	BRCA1_NonTargetingControlGuideForHuman0412	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
1	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	0	SETD2	SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352_SETD2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
2	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	BRCA1	BRCA1_chr17_41256141	0	NonTargetingControlGuideForHuman0362	BRCA1_NonTargetingControlGuideForHuman0362	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	153	321
3	NonTargetingControlGuideForHuman0352__NonTargetingControlGuideForHuman0412	0	NonTargetingControlGuideForHuman0352	0	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0352__NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0352__NonTargetingControlGuideForHuman0412	321	377	370	160
4	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	PIK3R1	PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412__PIK3R1	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
5	ALK_chr2_30143522__ALK_chr6_4415238	ALK	ALK_chr2_30143522	ALK	ALK_chr6_4415238	ALK__ALK	ALK_chr2_30143522__ALK_chr6_4415238	1536	714	2916	1391
"""

        expected_output_str = """construct_id	target_a_id	probe_a_id	target_b_id	probe_b_id	target_pair_id	probe_pair_id	A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
0	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	0	NonTargetingControlGuideForHuman0412	BRCA1_NonTargetingControlGuideForHuman0412	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
1	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	0	SETD2	SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352_SETD2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
2	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	BRCA1	BRCA1_chr17_41256141	0	NonTargetingControlGuideForHuman0362	BRCA1_NonTargetingControlGuideForHuman0362	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	153	321
4	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	PIK3R1	PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412__PIK3R1	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""

        scoring_data = TestScoringData.help_make_scoring_data(input_str)
        scoring_data._remove_constructs_w_both_probes_targeting_same_target_id()

        expected_output = ns_help.help_make_df(expected_output_str)
        self.assertTrue(expected_output.equals(scoring_data._input_counts_and_annotation_df))

    def test__get_data_only_df(self):
        input_str = """construct_id	target_a_id	probe_a_id	target_b_id	probe_b_id	target_pair_id	probe_pair_id	A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
0	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	0	NonTargetingControlGuideForHuman0412	BRCA1_NonTargetingControlGuideForHuman0412	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
1	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	0	SETD2	SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352_SETD2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
2	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	BRCA1	BRCA1_chr17_41256141	0	NonTargetingControlGuideForHuman0362	BRCA1_NonTargetingControlGuideForHuman0362	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	153	321
4	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	PIK3R1	PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412__PIK3R1	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""

        expected_output_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	153	321
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""
        scoring_data = TestScoringData.help_make_scoring_data(input_str)
        real_output_df = scoring_data._get_data_only_df()

        expected_output = ns_help.help_make_df(expected_output_str)
        self.assertTrue(expected_output.equals(real_output_df))

    def test__add_pseudocount(self):
        input_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	0	321
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""

        expected_output_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	1	48	1	1
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	1	321
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""
        data_only_df = pandas.read_csv(io.StringIO(input_str), sep="\t")
        real_output_df = ns_test.ScoringData._add_pseudocount(data_only_df, 1)

        expected_output = ns_help.help_make_df(expected_output_str)
        self.assertTrue(expected_output.equals(real_output_df))

    def test__add_pseudocount_custom_value(self):
        input_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	0	321
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""

        expected_output_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	2	48	2	2
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	2	321
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""
        data_only_df = pandas.read_csv(io.StringIO(input_str), sep="\t")
        real_output_df = ns_test.ScoringData._add_pseudocount(data_only_df, 2)

        expected_output = ns_help.help_make_df(expected_output_str)
        self.assertTrue(expected_output.equals(real_output_df))

    def test__calc_total_counts_per_sample(self):
        input_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	1	48	1	1
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	1	321
"""
        expected_output_keys = ["A549CV4_T21_1",
                                "A549CV4_T21_2",
                                "A549CV4_T28_1",
                                "A549CV4_T28_2"]
        expected_output_vals = [1375, 1168, 557, 1177]

        data_only_df = pandas.read_csv(io.StringIO(input_str), sep="\t")
        real_output_series = ns_test.ScoringData._calc_total_counts_per_sample(data_only_df)

        expected_output = pandas.Series(data=expected_output_vals, index=expected_output_keys)
        self.assertTrue(expected_output.equals(real_output_series))

    def test__calc_log2_fractions(self):
        input_pseudocounted_data_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	1	48	1	1
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	153	321
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	244	1341	100	2356
"""

        input_abundance_keys = ["A549CV4_T21_1",
                                "A549CV4_T21_2",
                                "A549CV4_T28_1",
                                "A549CV4_T28_2"]
        input_abundance_vals = [8292426, 8347478, 9995195, 11473098]

        expected_output_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	-22.983363	-17.407946	-23.252803	-23.451752
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	-12.777570	-13.549965	-14.136459	-13.711971
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	-15.390906	-14.264989	-15.995415	-15.125322
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	-15.052625	-12.603815	-16.608947	-12.249628
"""
        data_only_df = pandas.read_csv(io.StringIO(input_pseudocounted_data_str), sep="\t")
        abundances_series = pandas.Series(data=input_abundance_vals, index=input_abundance_keys)
        real_output_df = ns_test.ScoringData._calc_log2_fractions(data_only_df, abundances_series)

        expected_output = ns_help.help_make_df(expected_output_str)
        self.assertTrue(expected_output.round(decimals=6).equals(real_output_df.round(decimals=6)))

    def test__convert_abundance_thresholds_units(self):
        input_str = """sampleName	log2CountsThresh
A549CV4_T21_1	4.525
A549CV4_T21_2	4.325000000000001
A549CV4_T28_1	3.9250000000000003
A549CV4_T28_2	3.4750000000000005
"""
        # Note index_col param: first col should be used as index of resulting dataframe
        log2_counts_thresholds_per_sample_df = pandas.read_csv(io.StringIO(input_str), sep="\t", index_col=0)

        input_abundance_keys = ["A549CV4_T21_1",
                                "A549CV4_T21_2",
                                "A549CV4_T28_1",
                                "A549CV4_T28_2"]
        input_abundance_vals = [8292426, 8347478, 9995195, 11473098]
        total_counts_per_sample_series = pandas.Series(data=input_abundance_vals, index=input_abundance_keys)

        expected_output_keys = ["A549CV4_T21_1",
                                "A549CV4_T21_2",
                                "A549CV4_T28_1",
                                "A549CV4_T28_2"]
        expected_output_vals = [-18.45836, -18.66791, -19.32780, -19.97675]

        real_output_series = ns_test.ScoringData._convert_abundance_thresholds_units(
            total_counts_per_sample_series, log2_counts_thresholds_per_sample_df)

        expected_output_series = pandas.Series(data=expected_output_vals, index=expected_output_keys)
        self.assertTrue(expected_output_series.round(decimals=5).equals(real_output_series.round(decimals=5)))

    def test__convert_abundance_thresholds_units_error(self):
        input_str = """sampleName	log2CountsThresh
A549CV4_T21_1	4.525
A549CV4_T21_3	4.325000000000001
A549CV4_T28_1	3.9250000000000003
A549CV4_T28_3	3.4750000000000005
"""
        # Note index_col param: first col should be used as index of resulting dataframe
        expected_log2_counts_thresholds_per_sample_df = pandas.read_csv(io.StringIO(input_str), sep="\t", index_col=0)

        input_abundance_keys = ["A549CV4_T21_1",
                                "A549CV4_T21_2",
                                "A549CV4_T28_1",
                                "A549CV4_T28_2"]
        input_abundance_vals = [8292426, 8347478, 9995195, 11473098]
        total_counts_per_sample_series = pandas.Series(data=input_abundance_vals, index=input_abundance_keys)

        with self.assertRaises(ValueError):
            ns_test.ScoringData._convert_abundance_thresholds_units(
                total_counts_per_sample_series, expected_log2_counts_thresholds_per_sample_df)

    def test__get_column_headers_for_replicate(self):
        input_str = """construct_id	target_a_id	probe_a_id	target_b_id	probe_b_id	target_pair_id	probe_pair_id	A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
0	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	0	NonTargetingControlGuideForHuman0412	BRCA1_NonTargetingControlGuideForHuman0412	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
1	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	0	SETD2	SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352_SETD2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
2	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	BRCA1	BRCA1_chr17_41256141	0	NonTargetingControlGuideForHuman0362	BRCA1_NonTargetingControlGuideForHuman0362	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	153	321
4	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	PIK3R1	PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412__PIK3R1	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""
        expected_output = ["A549CV4_T21_2", "A549CV4_T28_2"]
        scoring_data = TestScoringData.help_make_scoring_data(input_str)
        real_output = scoring_data._get_column_headers_for_replicate(2)
        self.assertListEqual(expected_output, real_output)

    def test__extract_replicate_data_list(self):
        input_str = """construct_id	target_a_id	probe_a_id	target_b_id	probe_b_id	target_pair_id	probe_pair_id	A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
0	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	0	NonTargetingControlGuideForHuman0412	BRCA1_NonTargetingControlGuideForHuman0412	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
1	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	0	SETD2	SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352_SETD2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
2	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	BRCA1	BRCA1_chr17_41256141	0	NonTargetingControlGuideForHuman0362	BRCA1_NonTargetingControlGuideForHuman0362	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	153	321
4	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	PIK3R1	PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412__PIK3R1	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""

        input_log2_fractions_str = """A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	-5.67242534	-6.09759324	-10.09143539	-10.67065625
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	4.53336791	-2.23961225	-0.97509143	-0.93087564
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	1.9200317	-2.95463529	-2.83404754	-2.34422676
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	2.55157633	-3.4828834	-2.39099567	-0.97542796
"""
        input_log2_fractions_df = pandas.read_csv(io.StringIO(input_log2_fractions_str), sep="\t")

        input_log2_thresholds_keys = ["A549CV4_T21_1",
                                "A549CV4_T21_2",
                                "A549CV4_T28_1",
                                "A549CV4_T28_2"]
        input_log2_thresholds_vals = [-18.45836, -18.66791, -19.32780, -19.97675]
        input_log2_thresholds_series = pandas.Series(data=input_log2_thresholds_vals, index=input_log2_thresholds_keys)

        scoring_data = TestScoringData.help_make_scoring_data(input_str)
        scoring_data._timepoints_list = [21, 28]
        scoring_data._replicates_list = [1, 2]


        expected_rep1_log2_fractions_str = """A549CV4_T21_1	A549CV4_T28_1
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	-5.67242534	-10.09143539
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	4.53336791	-0.97509143
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	1.9200317	-2.83404754
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	2.55157633	-2.39099567
"""

        expected_rep2_log2_fractions_str = """A549CV4_T21_2	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	-6.09759324	-10.67065625
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	-2.23961225	-0.93087564
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	-2.95463529	-2.34422676
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	-3.4828834	-0.97542796
"""
        expected_output = [
            ns_per_rep.PerReplicateData(1, pandas.read_csv(io.StringIO(expected_rep1_log2_fractions_str), sep="\t"),
                                     pandas.Series(data=[-18.45836, -19.32780],
                                                   index=["A549CV4_T21_1","A549CV4_T28_1"])),
            ns_per_rep.PerReplicateData(2, pandas.read_csv(io.StringIO(expected_rep2_log2_fractions_str), sep="\t"),
                                     pandas.Series(data=[-18.66791, -19.97675],
                                                   index=["A549CV4_T21_2","A549CV4_T28_2"]))
        ]

        real_output = scoring_data._extract_replicate_data_list(input_log2_fractions_df,
                                                                    input_log2_thresholds_series)
        self.assertEqual(len(expected_output), len(real_output))
        for i in range(len(expected_output)):
            expected_item = expected_output[i]
            real_item = real_output[i]
            self.assertEqual(expected_item._replicate_num, real_item._replicate_num)
            self.assertTrue(expected_item._log2_fractions_thresholds_by_timepoints_series.round(decimals=6).equals(
                real_item._log2_fractions_thresholds_by_timepoints_series.round(decimals=6)))
            self.assertTrue(expected_item._log2_fractions_by_constructs_by_timepoints_df.round(decimals=6).equals(
                real_item._log2_fractions_by_constructs_by_timepoints_df.round(decimals=6)))

    def test_generate_per_replicate_data_list(self):
        input_str = """construct_id	target_a_id	probe_a_id	target_b_id	probe_b_id	target_pair_id	probe_pair_id	A549CV4_T21_1	A549CV4_T21_2	A549CV4_T28_1	A549CV4_T28_2
0	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	0	NonTargetingControlGuideForHuman0412	BRCA1_NonTargetingControlGuideForHuman0412	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	0	48	0	0
1	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	0	SETD2	SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352_SETD2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	1181	696	555	855
2	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	BRCA1	BRCA1_chr17_41256141	0	NonTargetingControlGuideForHuman0362	BRCA1_NonTargetingControlGuideForHuman0362	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	193	424	153	321
4	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	PIK3R1	PIK3R1_chr5_67576410	NonTargetingControlGuideForHuman0412__PIK3R1	NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	299	294	208	829
"""

        input_log2_counts_thresholds_per_sample_str = """sampleName	log2CountsThresh
A549CV4_T21_1	4.525
A549CV4_T21_2	4.325000000000001
A549CV4_T28_1	3.9250000000000003
A549CV4_T28_2	3.4750000000000005
"""
        # Note index_col param: first col should be used as index of resulting dataframe
        log2_counts_thresholds_per_sample_df = pandas.read_csv(io.StringIO(input_log2_counts_thresholds_per_sample_str),
                                                               sep="\t", index_col=0)

        scoring_data = TestScoringData.help_make_scoring_data(input_str)
        scoring_data._timepoints_list = [21, 28]
        scoring_data._replicates_list = [1, 2]

        expected_rep1_log2_fractions_str = """A549CV4_T21_1	A549CV4_T28_1
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	-10.709084	-9.840778
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	-0.503291	-0.724434
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	-3.116627	-2.58339
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	-2.485082	-2.140338
"""

        expected_rep2_log2_fractions_str = """A549CV4_T21_2	A549CV4_T28_2
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	-4.928765	-10.970106
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	-1.070784	-1.230325
BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	-1.785807	-2.643676
NonTargetingControlGuideForHuman0412__PIK3R1_chr5_67576410	-2.314055	-1.274878
"""

        expected_output = [
            ns_per_rep.PerReplicateData(1, pandas.read_csv(io.StringIO(expected_rep1_log2_fractions_str), sep="\t"),
                                     pandas.Series(data=[-6.184084, -5.915778],
                                                   index=pandas.Index(data=["A549CV4_T21_1", "A549CV4_T28_1"],
                                                                      name="sampleName"))),
            ns_per_rep.PerReplicateData(2, pandas.read_csv(io.StringIO(expected_rep2_log2_fractions_str), sep="\t"),
                                     pandas.Series(data=[-6.188728, -7.495106],
                                                   index=pandas.Index(data=["A549CV4_T21_2", "A549CV4_T28_2"],
                                                                      name="sampleName")))
        ]

        real_output = scoring_data.generate_per_replicate_data_list(log2_counts_thresholds_per_sample_df)

        self.assertEqual(len(expected_output), len(real_output))
        for i in range(len(expected_output)):
            expected_item = expected_output[i]
            real_item = real_output[i]
            self.assertEqual(expected_item._replicate_num, real_item._replicate_num)
            self.assertTrue(expected_item._log2_fractions_thresholds_by_timepoints_series.round(decimals=6).equals(
                real_item._log2_fractions_thresholds_by_timepoints_series.round(decimals=6)))
            self.assertTrue(expected_item._log2_fractions_by_constructs_by_timepoints_df.round(decimals=6).equals(
                real_item._log2_fractions_by_constructs_by_timepoints_df.round(decimals=6)))

    def test_main_functionality(self):
        input_data_df = ns_help.access_test_file_as_df("TestSet8_timepoint_counts_new.txt", col_index_to_use_as_index=None)
        input_abundance_df = ns_help.access_test_file_as_df("TestSet8_abundance_thresholds.txt")

        expected_output = [ns_help.TestPerReplicateData.help_make_per_rep_data(get_rep_1=True),
                           ns_help.TestPerReplicateData.help_make_per_rep_data(get_rep_1=False)]

        scoring_data = ns_test.ScoringData(input_data_df, "NonTargeting")
        real_output = scoring_data.generate_per_replicate_data_list(input_abundance_df)

        self.assertEqual(len(expected_output), len(real_output))
        for i in range(len(expected_output)):
            expected_item = expected_output[i]
            real_item = real_output[i]

            self.assertEqual(expected_item._replicate_num, real_item._replicate_num)
            self.assertEqual(expected_item._log2_fractions_thresholds_by_timepoints_series.round(decimals=5).to_csv(None),
                            real_item._log2_fractions_thresholds_by_timepoints_series.round(decimals=5).to_csv(None))
            # TODO: deal with creating comparable expected log2 fractions df, including index names (note 0 prefix in
            # Roman's code, not in mine) and associated sort order, as well as column names
            self.assertEqual(expected_item._log2_fractions_by_constructs_by_timepoints_df.round(decimals=5).to_csv(None),
                            real_item._log2_fractions_by_constructs_by_timepoints_df.round(decimals=5).to_csv(None))
