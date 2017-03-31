# standard libraries
import io
import unittest

# project-specific libraries
import ccbb_mali_dual_crispr.dual_crispr.construct_file_extracter as ns_extracter

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
    # region extract_construct_and_grna_info
    def test_extract_construct_and_grna_info(self):
        construct_input = io.StringIO("""# library_name = CV4
# min_trimmed_grna_len = 19
# max_trimmed_grna_len = 21
construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq	FinalSequence	Gene_A_chr	Gene_A_pos	Gene_B_chr	Gene_B_pos
NonTargetingControlGuideForHuman0412__FLT3_chr13_28636131	NonTargetingControlGuideForHuman0412_NA	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	FLT3-5	FLT3_chr13_28636131	GCGAGGCGCGCCGCTCCAGG	NonTargetingControlGuideForHuman0412_NA__FLT3-5	tatatatcttgtggaaaggacgaaacACCGGCACGCTGTACAGACGACAAGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCGAGGCGCGCCGCTCCAGGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
HDAC1_chr1_32757816__NonTargetingControlGuideForHuman0412	HDAC1-1	HDAC1_chr1_32757816	gcgctcgcgcccggacgcgg	NonTargetingControlGuideForHuman0412_NA	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	HDAC1-1__NonTargetingControlGuideForHuman0412_NA	tatatatcttgtggaaaggacgaaacACCGGACCGACTGACGGTAGGGACGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
FGFR2_chr10_123298215__NonTargetingControlGuideForHuman0412	FGFR2-13	FGFR2_chr10_123298215	gcgcggccgccACAAAGCTC	NonTargetingControlGuideForHuman0412_NA	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	FGFR2-13__NonTargetingControlGuideForHuman0412_NA	tatatatcttgtggaaaggacgaaacACCGgcgcggccgccACAAAGCTCGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
""")
        # note that: list items are in sorted order by probe name, NonTargetingControlGuideForHuman0412 only occurs once
        # although it occurs multiple times in input, and the two seqs that weren't already all upper case have been
        # uppercased.
        expected_grna_names_and_seqs = [("FGFR2_chr10_123298215", "GCGCGGCCGCCACAAAGCTC"),  # note: converted to upper
                                        ("FLT3_chr13_28636131", "GCGAGGCGCGCCGCTCCAGG"),
                                        ("HDAC1_chr1_32757816", "GCGCTCGCGCCCGGACGCGG"),  # note: converted to upper
                                        ("NonTargetingControlGuideForHuman0412", "GCACGCTGTACAGACGACAA")]

        expected_construct_names = ["NonTargetingControlGuideForHuman0412__FLT3_chr13_28636131",
                                    "HDAC1_chr1_32757816__NonTargetingControlGuideForHuman0412",
                                    "FGFR2_chr10_123298215__NonTargetingControlGuideForHuman0412"]

        construct_names, output_grna_names_and_seqs = ns_extracter.extract_construct_and_grna_info(construct_input)
        self.assertListEqual(expected_grna_names_and_seqs, output_grna_names_and_seqs)
        self.assertListEqual(expected_construct_names, construct_names)

    # endregion

    def test__read_in_construct_table(self):
        construct_input = io.StringIO("""# library_name = CV4
# min_trimmed_grna_len = 19
# max_trimmed_grna_len = 21
construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq	FinalSequence
NonTargetingControlGuideForHuman0412__FLT3_chr13_28636131	NonTargetingControlGuideForHuman0412_NA	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	FLT3-5	FLT3_chr13_28636131	GCGAGGCGCGCCGCTCCAGG	NonTargetingControlGuideForHuman0412_NA__FLT3-5	tatatatcttgtggaaaggacgaaacACCGGCACGCTGTACAGACGACAAGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCGAGGCGCGCCGCTCCAGGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
HDAC1_chr1_32757816__NonTargetingControlGuideForHuman0412	HDAC1-1	HDAC1_chr1_32757816	gcgctcgcgcccggacgcgg	NonTargetingControlGuideForHuman0412_NA	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	HDAC1-1__NonTargetingControlGuideForHuman0412_NA	tatatatcttgtggaaaggacgaaacACCGGACCGACTGACGGTAGGGACGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
FGFR2_chr10_123298215__NonTargetingControlGuideForHuman0412	FGFR2-13	FGFR2_chr10_123298215	gcgcggccgccACAAAGCTC	NonTargetingControlGuideForHuman0412_NA	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	FGFR2-13__NonTargetingControlGuideForHuman0412_NA	tatatatcttgtggaaaggacgaaacACCGgcgcggccgccACAAAGCTCGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
""")

        expected_col_names = ["construct_id", "target_a_id","probe_a_id","probe_a_seq","target_b_id","probe_b_id",
                              "probe_b_seq",7, 8, 9]

        real_output = ns_extracter._read_in_construct_table(construct_input)
        self.assertListEqual(expected_col_names, list(real_output.columns.values))

