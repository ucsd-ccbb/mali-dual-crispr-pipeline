# standard libraries
import io
import unittest

# third-party libraries
import pandas

# project-specific libraries
from ccbbucsd.malicrispr.construct_file_extracter import get_target_id_header, \
    get_probe_id_header, get_construct_header, get_target_pair_id_header, \
    get_probe_pair_id_header

from ccbbucsd.malicrispr.scoring_prep import  \
    _validate_and_standardize_timepoint, _validate_expt_structure, \
    _validate_and_decompose_count_header, _validate_and_parse_data_column_headers, \
    _generate_scoring_friendly_annotation, _validate_and_standardize_replicate, \
    _validate_and_recompose_count_header

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
    # region _validate_and_decompose_count_header
    def test__validate_and_decompose_count_header_valid(self):
        header = "PGP1MV4_t21_2"
        expected_output = ("PGP1MV4", 21, 2)
        real_output = _validate_and_decompose_count_header(header)
        self.assertEqual(expected_output, real_output)

        real_output2 = _validate_and_decompose_count_header("PGP1-MV4_t21_2_S9_trimmed53_len_filtered_counts")
        self.assertEqual(("PGP1-MV4", 21, 2), real_output2)

        real_output3 = _validate_and_decompose_count_header("PGP1-MV4_t21_2")
        self.assertEqual(("PGP1-MV4", 21, 2), real_output2)

    def test__validate_and_decompose_count_header_invalid(self):
        with self.assertRaises(ValueError):
            _validate_and_decompose_count_header("PGP1_MV4_t21_2_S9_trimmed53_len_filtered_counts")

        with self.assertRaises(ValueError):
            _validate_and_decompose_count_header("PGP1_MV4_t21_2")

        with self.assertRaises(ValueError):
            _validate_and_decompose_count_header("PGPrep1")

    # end region

    # region _validate_and_standardize_timepoint
    def test__validate_and_standardize_timepoint_valid(self):
        expected_output = 40
        real_output = _validate_and_standardize_timepoint("t40")
        self.assertEqual(expected_output, real_output)

        real_output_2 = _validate_and_standardize_timepoint("T40")
        self.assertEqual(expected_output, real_output_2)

    def test__validate_and_standardize_timepoint_invalid_letter(self):
        with self.assertRaises(ValueError):
            # any letter other than "t" or "T" at beginnning of
            # timepoint should fail
            _validate_and_standardize_timepoint("d40")

    def test__validate_and_standardize_timepoint_invalid_number(self):
        with self.assertRaises(ValueError):
            _validate_and_standardize_timepoint("test40")

        with self.assertRaises(ValueError):
            _validate_and_standardize_timepoint("t-40")

        with self.assertRaises(ValueError):
            _validate_and_standardize_timepoint("t4.1")

    # end region

    # region _validate_and_standardize_replicate
    def test__validate_and_standardize_replicate_digit(self):
        real_output = _validate_and_standardize_replicate("10")
        self.assertEqual(10, real_output)

    def test__validate_and_standardize_replicate_non_digit(self):
        input = "A"
        real_output = _validate_and_standardize_replicate(input)
        self.assertEqual(input, real_output)

        input_2 = "4.1"
        real_output_2 = _validate_and_standardize_replicate(input_2)
        self.assertEqual(input_2, real_output_2)

    # end region

    # region _validate_expt_structure
    def test__validate_expt_structure_valid(self):
        input = {
            "expt1": {
                "T0": {"1", "2", "3"},
                "T1": {"1", "2", "3"}}
            }

        _validate_expt_structure(input)

        # If we got this far, then the test
        # passed by definition
        self.assertTrue(True, True)

    def test__validate_expt_structure_diff_reps_for_timept_in_expt(self):
        input = {
            "expt1": {
                "T0": {"1", "2", "3"},
                "T1": {"1", "2"}}
        }

        with self.assertRaises(ValueError):
            _validate_expt_structure(input)

    def test__validate_expt_structure_multiple_expts(self):
        input = {
            "expt1": {
                "T0": {"1", "2", "3"},
                "T1": {"1", "2", "3"}},
            "expt2": {
                "T0": {"1", "2", "3"},
                "T1": {"1", "2", "3"}}
            }

        with self.assertRaises(ValueError):
            _validate_expt_structure(input)

    def test__validate_expt_structure_no_timepts(self):
        input = {
            "expt1": {}
        }

        with self.assertRaises(ValueError):
            _validate_expt_structure(input)

    def test__validate_expt_structure_no_reps(self):
        input = {
            "expt1": {
                "T0": {}},
        }

        with self.assertRaises(ValueError):
            _validate_expt_structure(input)

    # def test__validate_expt_structure_valid(self):
    #     input = {
    #         "expt1": {
    #             "T0": {"1", "2", "3"},
    #             "T1": {"1", "2", "3"}},
    #         "expt2": {
    #             "T0": {"1", "2", "3"},
    #             "T1": {"1", "2", "3"}}
    #         }
    #
    #     _validate_expt_structure(input)
    #
    #     # If we got this far, then the test
    #     # passed by definition
    #     self.assertTrue(True, True)
    #
    # def test__validate_expt_structure_diff_reps_for_timept_in_expt(self):
    #     input = {
    #         "expt1": {
    #             "T0": {"1", "2", "3"},
    #             "T1": {"1", "2"}},
    #         "expt2": {
    #             "T0": {"1", "2", "3"},
    #             "T1": {"1", "2"}}
    #     }
    #
    #     with self.assertRaises(ValueError):
    #         _validate_expt_structure(input)
    #
    # def test__validate_expt_structure_diff_timept_reps_for_expts(self):
    #     input = {
    #         "expt1": {
    #             "T0": {"1", "2", "3"},
    #             "T1": {"1", "2", "3"}},
    #         "expt2": {
    #             "T0": {"1", "2", "3"},
    #             "T1": {"1", "2", "3"},
    #             "T2": {"1", "2", "3"}},
    #     }
    #
    #     with self.assertRaises(ValueError):
    #         _validate_expt_structure(input)

    # end region

    # region _validate_and_parse_data_column_headers
    def test__validate_and_parse_data_column_headers_valid_num_reps(self):
        data_headers = ["A549-MV4_t28_1_S7_trimmed53_len_filtered_counts",
                        "A549-MV4_t14_2_S4_trimmed53_len_filtered_counts",
                        "A549-MV4_t28_2_S8_trimmed53_len_filtered_counts",
                        "A549-MV4_t3_2_S2_trimmed53_len_filtered_counts",
                        "A549-MV4_t20_2_S6_trimmed53_len_filtered_counts",
                        "A549-MV4_t14_1_S3_trimmed53_len_filtered_counts",
                        "A549-MV4_t3_1_S1_trimmed53_len_filtered_counts",
                        "A549-MV4_t20_1_S5_trimmed53_len_filtered_counts"]

        expected_output = [('A549MV4', 28, 1), ('A549MV4', 14, 2), ('A549MV4', 28, 2), ('A549MV4', 3, 2),
                           ('A549MV4', 20, 2), ('A549MV4', 14, 1), ('A549MV4', 3, 1), ('A549MV4', 20, 1)]

        real_output = _validate_and_parse_data_column_headers(data_headers, "A549MV4")
        self.assertEqual(expected_output, real_output)

    def test__validate_and_parse_data_column_headers_valid_nonnum_reps(self):
        data_headers = ["A549-MV4_t28_a_S7_trimmed53_len_filtered_counts",
                        "A549-MV4_t14_b_S4_trimmed53_len_filtered_counts",
                        "A549-MV4_t28_b_S8_trimmed53_len_filtered_counts",
                        "A549-MV4_t3_b_S2_trimmed53_len_filtered_counts",
                        "A549-MV4_t20_b_S6_trimmed53_len_filtered_counts",
                        "A549-MV4_t14_a_S3_trimmed53_len_filtered_counts",
                        "A549-MV4_t3_a_S1_trimmed53_len_filtered_counts",
                        "A549-MV4_t20_a_S5_trimmed53_len_filtered_counts"]

        expected_output = [('A549MV4', 28, 'a'), ('A549MV4', 14, 'b'), ('A549MV4', 28, 'b'), ('A549MV4', 3, 'b'),
                           ('A549MV4', 20, 'b'), ('A549MV4', 14, 'a'), ('A549MV4', 3, 'a'), ('A549MV4', 20, 'a')]

        real_output = _validate_and_parse_data_column_headers(data_headers, "A549MV4")
        self.assertEqual(expected_output, real_output)

    # end region

    # region _generate_scoring_friendly_annotation
    def test__generate_scoring_friendly_annotation(self):
        input_str = """0	construct_id	2	target_a_id	4	5	probe_a_id	probe_a_seq	target_b_id	9	10	probe_b_id	probe_b_seq
1	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	tatatatcttgtggaaaggacgaaacACCGTCTTGTGCTGACTTACCAGAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	BRCA1	chr17	41276018	BRCA1_chr17_41276018	TCTTGTGCTGACTTACCAGA	NonTargetingControlGuideForHuman0412			NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA
2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	tatatatcttgtggaaaggacgaaacACCGGCCATTCTAGTCCCGGCATAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGAACGTTAACTCTGAGCCTGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	NonTargetingControlGuideForHuman0352			NonTargetingControlGuideForHuman0352	GCCATTCTAGTCCCGGCATA	SETD2	chr3	47142972	SETD2_chr3_47142972	AACGTTAACTCTGAGCCTGA
3	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	tatatatcttgtggaaaggacgaaacACCGTGAACCCGAAAATCCTTCCTGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGAGTGATGCTTAGACTCCGTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	BRCA1	chr17	41256141	BRCA1_chr17_41256141	TGAACCCGAAAATCCTTCCT	NonTargetingControlGuideForHuman0362			NonTargetingControlGuideForHuman0362	GAGTGATGCTTAGACTCCGT"""
        expected_output_str = """	construct_id	target_a_id	probe_a_id	target_b_id	probe_b_id	target_pair_id	probe_pair_id
0	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	BRCA1_NonTargetingControlGuideForHuman0412	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412
1	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	NonTargetingControlGuideForHuman0352	SETD2	SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352_SETD2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972
2	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	BRCA1	BRCA1_chr17_41256141	NonTargetingControlGuideForHuman0362	NonTargetingControlGuideForHuman0362	BRCA1_NonTargetingControlGuideForHuman0362	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362
"""

        input_df = pandas.read_csv(io.StringIO(input_str), sep="\t")
        expected_output_df = pandas.read_csv(io.StringIO(expected_output_str), sep="\t")
        real_output_df = _generate_scoring_friendly_annotation(input_df)

        expected_column_headers = [get_construct_header(),
            get_probe_id_header("a"), get_probe_id_header("b"),
            get_target_id_header("a"), get_target_id_header("b")]
        self.assertEqual(expected_column_headers, real_output_df.columns.values.tolist())

        # these dfs aren't equal, but some of their columns should be:
        self.assertEqual(expected_output_df[get_construct_header()].tolist(),
                         real_output_df[get_construct_header()].tolist())
        self.assertEqual(expected_output_df[get_target_id_header("a")].tolist(),
                         real_output_df[get_target_id_header("a")].tolist())
        self.assertEqual(expected_output_df[get_target_id_header("b")].tolist(),
                         real_output_df[get_target_id_header("b")].tolist())
        self.assertEqual(expected_output_df[get_probe_id_header("a")].tolist(),
                         real_output_df[get_probe_id_header("a")].tolist())
        self.assertEqual(expected_output_df[get_probe_id_header("b")].tolist(),
                         real_output_df[get_probe_id_header("b")].tolist())

    # This test represents what I'm trying to refactor the scoring data prep code to accept; not there yet.
    #     def test__generate_scoring_friendly_annotation(self):
    #         input_str = """0	construct_id	2	target_a_id	4	5	probe_a_id	probe_a_seq	target_b_id	9	10	probe_b_id	probe_b_seq
    # 1	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	tatatatcttgtggaaaggacgaaacACCGTCTTGTGCTGACTTACCAGAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	BRCA1	chr17	41276018	BRCA1_chr17_41276018	TCTTGTGCTGACTTACCAGA	NonTargetingControlGuideForHuman0412			NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA
    # 2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	tatatatcttgtggaaaggacgaaacACCGGCCATTCTAGTCCCGGCATAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGAACGTTAACTCTGAGCCTGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	NonTargetingControlGuideForHuman0352			NonTargetingControlGuideForHuman0352	GCCATTCTAGTCCCGGCATA	SETD2	chr3	47142972	SETD2_chr3_47142972	AACGTTAACTCTGAGCCTGA
    # 3	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	tatatatcttgtggaaaggacgaaacACCGTGAACCCGAAAATCCTTCCTGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGAGTGATGCTTAGACTCCGTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	BRCA1	chr17	41256141	BRCA1_chr17_41256141	TGAACCCGAAAATCCTTCCT	NonTargetingControlGuideForHuman0362			NonTargetingControlGuideForHuman0362	GAGTGATGCTTAGACTCCGT"""
    #         expected_output_str = """	construct_id	target_a_id	probe_a_id	target_b_id	probe_b_id	target_pair_id	probe_pair_id
    # 0	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	BRCA1_NonTargetingControlGuideForHuman0412	BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412
    # 1	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	NonTargetingControlGuideForHuman0352	SETD2	SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352_SETD2	NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972
    # 2	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362	BRCA1	BRCA1_chr17_41256141	NonTargetingControlGuideForHuman0362	NonTargetingControlGuideForHuman0362	BRCA1_NonTargetingControlGuideForHuman0362	BRCA1_chr17_41256141__NonTargetingControlGuideForHuman0362
    # """
    #
    #         input_df = pandas.read_csv(io.StringIO(input_str), sep="\t")
    #         expected_output_df = pandas.read_csv(io.StringIO(expected_output_str), sep="\t")
    #         real_output_df = _generate_scoring_friendly_annotation(input_df)
    #
    #         expected_column_headers = [get_construct_header(), get_target_id_header("a"),
    #                                    get_probe_id_header("a"), get_target_id_header("b"), get_probe_id_header("b"),
    #                                    get_target_pair_id_header(), get_probe_pair_id_header()]
    #         self.assertEqual(expected_column_headers, real_output_df.columns.values.tolist())
    #
    #         # these dfs aren't equal, but some of their columns should be:
    #         self.assertEqual(expected_output_df[get_construct_header()].tolist(),
    #                          real_output_df[get_construct_header()].tolist())
    #         self.assertEqual(expected_output_df[get_target_id_header("a")].tolist(),
    #                          real_output_df[get_target_id_header("a")].tolist())
    #         self.assertEqual(expected_output_df[get_target_id_header("b")].tolist(),
    #                          real_output_df[get_target_id_header("b")].tolist())
    #         self.assertEqual(expected_output_df[get_probe_id_header("a")].tolist(),
    #                          real_output_df[get_probe_id_header("a")].tolist())
    #         self.assertEqual(expected_output_df[get_probe_id_header("b")].tolist(),
    #                          real_output_df[get_probe_id_header("b")].tolist())
    #         self.assertEqual(expected_output_df[get_target_pair_id_header()].tolist(),
    #                          real_output_df[get_target_pair_id_header()].tolist())
    #         self.assertEqual(expected_output_df[get_probe_pair_id_header()].tolist(),
    #                          real_output_df[get_probe_pair_id_header()].tolist())

    # end region

    # region _validate_and_recompose_count_header
    def test__validate_and_recompose_count_header(self):
        input = ("A549MV4", 3, 1)
        expected_output = "A549MV4_T3_1"
        real_output = _validate_and_recompose_count_header(input)
        self.assertEqual(expected_output, real_output)

    def test__validate_and_recompose_count_header_error(self):
        input = ("A549MV4", 3)
        with self.assertRaises(ValueError):
            _validate_and_recompose_count_header(input)

    # end region