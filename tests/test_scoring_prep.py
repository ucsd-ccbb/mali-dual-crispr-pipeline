# standard libraries
import io
import unittest

# third-party libraries
import pandas

# project-specific libraries
import ccbb_mali_dual_crispr.dual_crispr.construct_file_extracter as ns_extractor
import ccbb_mali_dual_crispr.dual_crispr.scoring_prep as ns_test

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
    # region read_timepoint_from_standardized_count_header

    def test_read_timepoint_from_standardized_count_header_valid(self):
        real_output = ns_test.read_timepoint_from_standardized_count_header("PGP1MV4_t4_1", ["T","D"])
        self.assertEqual(4, real_output)

        real_output2 = ns_test.read_timepoint_from_standardized_count_header("PGP1MV4_d40_1", ["T","D"])
        self.assertEqual(40, real_output2)

    def test_read_timepoint_from_standardized_count_header_invalid_timept(self):
        with self.assertRaises(ValueError):
            ns_test.read_timepoint_from_standardized_count_header("PGP1MV4_q4_1", ["T","D"])

        with self.assertRaises(ValueError):
            ns_test.read_timepoint_from_standardized_count_header("PGP1MV4_test40_1", ["T","D"])

        with self.assertRaises(ValueError):
            ns_test.read_timepoint_from_standardized_count_header("PGP1MV4_t-40_1", ["T","D"])

        with self.assertRaises(ValueError):
            ns_test.read_timepoint_from_standardized_count_header("PGP1MV4_t4.1_1", ["T","D"])

    def test_read_timepoint_from_standardized_count_header_invalid_header_too_few_pieces(self):
        with self.assertRaises(ValueError):
            ns_test.read_timepoint_from_standardized_count_header("PGP1MV4-t4.1-1", ["T","D"])

    def test_read_timepoint_from_standardized_count_header_invalid_header_too_many_pieces(self):
        with self.assertRaises(ValueError):
            ns_test.read_timepoint_from_standardized_count_header("PGP1_MV4_t4.1_1", ["T","D"])

    # end region

    # region merge_and_annotate_counts
    def test_merge_and_annotate_counts(self):
        constructs_file_str = """# library_name = CV4
# min_trimmed_grna_len = 19
# max_trimmed_grna_len = 21
construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq	FinalSequence	Gene_A_chr	Gene_A_pos	Gene_B_chr	Gene_B_pos
NonTargetingControlGuideForHuman0412__FLT3_chr13_28636131	NonTargetingControlGuideForHuman0412_NA	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	FLT3-5	FLT3_chr13_28636131	GCGAGGCGCGCCGCTCCAGG	NonTargetingControlGuideForHuman0412_NA__FLT3-5	tatatatcttgtggaaaggacgaaacACCGGCACGCTGTACAGACGACAAGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCGAGGCGCGCCGCTCCAGGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
HDAC1_chr1_32757816__NonTargetingControlGuideForHuman0412	HDAC1-1	HDAC1_chr1_32757816	gcgctcgcgcccggacgcgg	NonTargetingControlGuideForHuman0412_NA	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	HDAC1-1__NonTargetingControlGuideForHuman0412_NA	tatatatcttgtggaaaggacgaaacACCGGACCGACTGACGGTAGGGACGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
FGFR2_chr10_123298215__NonTargetingControlGuideForHuman0412	FGFR2-13	FGFR2_chr10_123298215	gcgcggccgccACAAAGCTC	NonTargetingControlGuideForHuman0412_NA	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	FGFR2-13__NonTargetingControlGuideForHuman0412_NA	tatatatcttgtggaaaggacgaaacACCGgcgcggccgccACAAAGCTCGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
"""
        count_file_str = """construct_id	PGP1_MV4_t3_2_S9_trimmed53_len_filtered_counts	PGP1-MV4_t21_2_S9_trimmed53_len_filtered_counts
NonTargetingControlGuideForHuman0412__FLT3_chr13_28636131	40	1874
HDAC1_chr1_32757816__NonTargetingControlGuideForHuman0412	0	37
FGFR2_chr10_123298215__NonTargetingControlGuideForHuman0412	1938	1736
"""
        merge_df_str = """construct_id	probe_a_id	probe_b_id	target_a_id	target_b_id	PGP1MV4_T3_2	PGP1MV4_T21_2
NonTargetingControlGuideForHuman0412__FLT3_chr13_28636131	FLT3_chr13_28636131	NonTargetingControlGuideForHuman0412	FLT3-5	NonTargetingControlGuideForHuman0412_NA	40	1874
HDAC1_chr1_32757816__NonTargetingControlGuideForHuman0412	HDAC1_chr1_32757816	NonTargetingControlGuideForHuman0412	HDAC1-1	NonTargetingControlGuideForHuman0412_NA	0	37
FGFR2_chr10_123298215__NonTargetingControlGuideForHuman0412	FGFR2_chr10_123298215	NonTargetingControlGuideForHuman0412	FGFR2-13	NonTargetingControlGuideForHuman0412_NA	1938	1736
"""
        merge_df_file_obj = io.StringIO(merge_df_str)
        expected_output = pandas.read_table(merge_df_file_obj)

        construct_file_obj = io.StringIO(constructs_file_str)
        count_file_obj = io.StringIO(count_file_str)
        real_output = ns_test.merge_and_annotate_counts([count_file_obj], construct_file_obj, "PGP1MV4", ["T", "D"])
        self.assertTrue(expected_output.equals(real_output))

    # end region

    # region _get_orig_count_headers
    def test__get_orig_count_headers(self):
        count_file_str = """construct_id\tPGP1_MV4_t1_2_S9_trimmed53_len_filtered_counts\tPGP1-MV4_t21_2_S9_trimmed53_len_filtered_counts
NonTargetingControlGuideForHuman0412__FLT3_chr13_28636131\t40\t1874
HDAC1_chr1_32757816__NonTargetingControlGuideForHuman0412\t0\t37
FGFR2_chr10_123298215__NonTargetingControlGuideForHuman0412\t1938\t1736
"""
        expected_output = ["PGP1_MV4_t1_2_S9_trimmed53_len_filtered_counts",
                           "PGP1-MV4_t21_2_S9_trimmed53_len_filtered_counts"]

        count_file_obj = io.StringIO(count_file_str)
        count_file_df = pandas.read_table(count_file_obj, sep="\t")
        real_output = ns_test._get_orig_count_headers(count_file_df)
        self.assertListEqual(expected_output, real_output)

    # end region

    # region _validate_and_standardize_header_pieces
    def test__validate_and_standardize_header_pieces_valid(self):
        header = "PGP1MV4_t21_2"
        expected_output = ("PGP1MV4", 21, 2)
        real_output = ns_test._validate_and_standardize_header_pieces(header, "PGP1MV4", ["T","D"])
        self.assertEqual(expected_output, real_output)

        real_output2 = ns_test._validate_and_standardize_header_pieces(
            "PGP1-MV4_t21_2_S9_trimmed53_len_filtered_counts", "PGP1MV4", ["T","D"])
        self.assertEqual(("PGP1MV4", 21, 2), real_output2)

        real_output3 = ns_test._validate_and_standardize_header_pieces("PGP1-MV4_t21_2", "PGP1MV4",
                                                                       ["T","D"])
        self.assertEqual(("PGP1MV4", 21, 2), real_output3)

        real_output4 = ns_test._validate_and_standardize_header_pieces("PGP1_MV4_t21_2", "PGP1MV4",
                                                                       ["T","D"])
        self.assertEqual(("PGP1MV4", 21, 2), real_output4)

    def test___validate_and_standardize_header_pieces_invalid(self):
        with self.assertRaises(ValueError):
            ns_test._validate_and_standardize_header_pieces("PGP1_MV4_t_21_2_S9_trimmed53_len_filtered_counts",
                                                            "PGP1MV4", ["T","D"])
        with self.assertRaises(ValueError):
            ns_test._validate_and_standardize_header_pieces("PGPrep1", "PGP1MV4", ["T","D"])

    # end region

    # region _validate_and_standardize_timept_and_replicate
    def test__validate_and_standardize_timept_and_replicate_valid(self):
        real_timept, real_replicate = ns_test._validate_and_standardize_timept_and_replicate("PGP1MV4_t4_1", ["T", "D"])
        self.assertEqual(4, real_timept)
        self.assertEqual(1, real_replicate)

    def test__validate_and_standardize_timept_and_replicate_invalid_too_few_pieces(self):
        with self.assertRaises(ValueError):
            ns_test._validate_and_standardize_timept_and_replicate("PGP1MV4-t4-1", ["T", "D"])

    # I am not testing all the functionality of _validate_and_standardize_timepoint and
    # _validate_and_standardize_replicate *again* here ... this is already a whitebox test (since I'm testing a
    # private function) so anyone who modifies it should know it references those other two.

    # end region

    # region _validate_and_standardize_timepoint
    def test__validate_and_standardize_timepoint_valid(self):
        expected_output = 40
        real_output = ns_test._validate_and_standardize_timepoint("d40", ["T","D"])
        self.assertEqual(expected_output, real_output)

        real_output_2 = ns_test._validate_and_standardize_timepoint("T40", ["T","D"])
        self.assertEqual(expected_output, real_output_2)

    def test__validate_and_standardize_timepoint_invalid_letter(self):
        with self.assertRaises(ValueError):
            ns_test._validate_and_standardize_timepoint("r40", ["T","D"])

    def test__validate_and_standardize_timepoint_invalid_number(self):
        with self.assertRaises(ValueError):
            ns_test._validate_and_standardize_timepoint("test40", ["T","D"])

        with self.assertRaises(ValueError):
            ns_test._validate_and_standardize_timepoint("t-40", ["T","D"])

        with self.assertRaises(ValueError):
            ns_test._validate_and_standardize_timepoint("t4.1", ["T","D"])

    # end region

    # region _validate_and_standardize_replicate
    def test__validate_and_standardize_replicate_digit(self):
        real_output = ns_test._validate_and_standardize_replicate("10")
        self.assertEqual(10, real_output)

    def test__validate_and_standardize_replicate_non_digit(self):
        input = "A"
        real_output = ns_test._validate_and_standardize_replicate(input)
        self.assertEqual(input, real_output)

        input_2 = "4.1"
        real_output_2 = ns_test._validate_and_standardize_replicate(input_2)
        self.assertEqual(input_2, real_output_2)

    # end region

    # region _validate_expt_structure
    def test__validate_expt_structure_valid(self):
        input = {
            "expt1": {
                "T0": {"1", "2", "3"},
                "T1": {"1", "2", "3"}}
            }

        ns_test._validate_expt_structure(input)

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
            ns_test._validate_expt_structure(input)

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
            ns_test._validate_expt_structure(input)

    def test__validate_expt_structure_no_timepts(self):
        input = {
            "expt1": {}
        }

        with self.assertRaises(ValueError):
            ns_test._validate_expt_structure(input)

    def test__validate_expt_structure_no_reps(self):
        input = {
            "expt1": {
                "T0": {}},
        }

        with self.assertRaises(ValueError):
            ns_test._validate_expt_structure(input)

    # end region

    # region _validate_and_standardize_count_headers
    def test__validate_and_standardize_count_headers_valid_num_reps(self):
        input_prefixes = ["T","D"]
        data_headers = ["A549-MV4_t28_1_S7_trimmed53_len_filtered_counts",
                        "A549-MV4_t14_2_S4_trimmed53_len_filtered_counts",
                        "A549-MV4_t28_2_S8_trimmed53_len_filtered_counts",
                        "A549-MV4_t3_2_S2_trimmed53_len_filtered_counts",
                        "A549-MV4_t20_2_S6_trimmed53_len_filtered_counts",
                        "A549-MV4_t14_1_S3_trimmed53_len_filtered_counts",
                        "A549-MV4_t3_1_S1_trimmed53_len_filtered_counts",
                        "A549-MV4_t20_1_S5_trimmed53_len_filtered_counts"]

        expected_output = ['A549MV4_T28_1', 'A549MV4_T14_2', 'A549MV4_T28_2', 'A549MV4_T3_2',
                           'A549MV4_T20_2', 'A549MV4_T14_1', 'A549MV4_T3_1', 'A549MV4_T20_1']

        real_output = ns_test._validate_and_standardize_count_headers(data_headers, "A549MV4", input_prefixes)
        self.assertEqual(expected_output, real_output)

    def test__validate_and_standardize_count_headers_valid_nonnum_reps(self):
        input_prefixes = ["T","D"]
        data_headers = ["A549-MV4_t28_a_S7_trimmed53_len_filtered_counts",
                        "A549-MV4_t14_b_S4_trimmed53_len_filtered_counts",
                        "A549-MV4_t28_b_S8_trimmed53_len_filtered_counts",
                        "A549-MV4_t3_b_S2_trimmed53_len_filtered_counts",
                        "A549-MV4_t20_b_S6_trimmed53_len_filtered_counts",
                        "A549-MV4_t14_a_S3_trimmed53_len_filtered_counts",
                        "A549-MV4_t3_a_S1_trimmed53_len_filtered_counts",
                        "A549-MV4_t20_a_S5_trimmed53_len_filtered_counts"]

        expected_output = ['A549MV4_T28_a', 'A549MV4_T14_b', 'A549MV4_T28_b', 'A549MV4_T3_b',
                           'A549MV4_T20_b', 'A549MV4_T14_a', 'A549MV4_T3_a', 'A549MV4_T20_a']

        real_output = ns_test._validate_and_standardize_count_headers(data_headers, "A549MV4", input_prefixes)
        self.assertEqual(expected_output, real_output)

    def test__validate_and_standardize_count_headers_invalid(self):
        input_prefixes = ["T","D"]
        data_headers = ["A549-MV4_t28_1_S7_trimmed53_len_filtered_counts",
                        "A549-MV4-t14-2_S4_trimmed53_len_filtered_counts",
                        "A549-MV4_q28_2_S8_trimmed53_len_filtered_counts",
                        "A549-MV4_t3_2_S2_trimmed53_len_filtered_counts",
                        "A549-MV4_t20_2_S6_trimmed53_len_filtered_counts",
                        "A549-MV4_t14.1_1_S3_trimmed53_len_filtered_counts",
                        "A549-MV4_t3_1_S1_trimmed53_len_filtered_counts",
                        "A549-MV4_t20_1_S5_trimmed53_len_filtered_counts",
                        "A549rerun-MV4_t3_1_S2_trimmed53_len_filtered_counts",
                        "A549rerun-MV4_t20_1_S3_trimmed53_len_filtered_counts"]

        expected_error_msg = """The following error(s) were detected during count file parsing:
Column header 'A549-MV4-t14-2_S4_trimmed53_len_filtered_counts' separates on the '_' delimiter into the following 1 piece(s) instead of the expected 2: 2.
Time point 'q28' does not start with upper or lower case versions of any of the expected prefixes T, D.
Time point value '14.1' is not recognizable as a positive integer.
The following pair of column headers both appear to represent the same timepoint and replicate: 'A549rerun-MV4_t3_1_S2_trimmed53_len_filtered_counts', 'A549-MV4_t3_1_S1_trimmed53_len_filtered_counts'.  Please modify the inputs to remove this ambiguity.
The following pair of column headers both appear to represent the same timepoint and replicate: 'A549rerun-MV4_t20_1_S3_trimmed53_len_filtered_counts', 'A549-MV4_t20_1_S5_trimmed53_len_filtered_counts'.  Please modify the inputs to remove this ambiguity."""
        with self.assertRaises(ValueError) as valerr:
            ns_test._validate_and_standardize_count_headers(data_headers, "A549MV4", input_prefixes)
        self.assertEqual(expected_error_msg, str(valerr.exception))

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
        real_output_df = ns_test._generate_scoring_friendly_annotation(input_df)

        expected_column_headers = [ns_extractor.get_construct_header(),
                                   ns_extractor.get_probe_id_header("a"), ns_extractor.get_probe_id_header("b"),
                                   ns_extractor.get_target_id_header("a"), ns_extractor.get_target_id_header("b")]
        self.assertEqual(expected_column_headers, real_output_df.columns.values.tolist())

        # these dfs aren't equal, but some of their columns should be:
        self.assertEqual(expected_output_df[ns_extractor.get_construct_header()].tolist(),
                         real_output_df[ns_extractor.get_construct_header()].tolist())
        self.assertEqual(expected_output_df[ns_extractor.get_target_id_header("a")].tolist(),
                         real_output_df[ns_extractor.get_target_id_header("a")].tolist())
        self.assertEqual(expected_output_df[ns_extractor.get_target_id_header("b")].tolist(),
                         real_output_df[ns_extractor.get_target_id_header("b")].tolist())
        self.assertEqual(expected_output_df[ns_extractor.get_probe_id_header("a")].tolist(),
                         real_output_df[ns_extractor.get_probe_id_header("a")].tolist())
        self.assertEqual(expected_output_df[ns_extractor.get_probe_id_header("b")].tolist(),
                         real_output_df[ns_extractor.get_probe_id_header("b")].tolist())

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

    # region _recompose_count_header

    def test__recompose_count_header(self):
        input_prefixes = ["T","t","D","d"]
        expected_output = "A549MV4_T3_1"
        real_output = ns_test._recompose_count_header("A549MV4", 3, 1, input_prefixes)
        self.assertEqual(expected_output, real_output)

    # end region

    def test__sort_headers(self):
        input_headers = ["PGP1MV4_T3_1", "PGP1MV4_T10_2", "PGP1MV4_T21_2", "PGP1MV4_T10_10"]
        expected_output = ["PGP1MV4_T3_1","PGP1MV4_T10_2","PGP1MV4_T10_10","PGP1MV4_T21_2"]
        real_output = ns_test._sort_headers(input_headers, "PGP1MV4", ["T","D"])
        self.assertEqual(expected_output, real_output)