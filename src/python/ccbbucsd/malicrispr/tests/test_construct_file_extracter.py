# standard libraries
import io
import unittest

# project-specific libraries
from ccbbucsd.malicrispr.construct_file_extracter import extract_construct_and_grna_info

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
    # region _extract_unique_probe_seq_pairs tests
    def test__extract_probes_from_construct_table(self):
        construct_input = io.StringIO("""	SequenceID	Gene_A	Gene_A_seq	Gene_B	Gene_B_seq	Grna_A	Grna_B	Final_id	Final_seq
1	NonTargetingControlGuideForHuman0412__FLT3_chr13_28636131	NonTargetingControlGuideForHuman0412_NA	GCACGCTGTACAGACGACAA	FLT3-5	GCGAGGCGCGCCGCTCCAGG	NonTargetingControlGuideForHuman0412	FLT3_chr13_28636131	NonTargetingControlGuideForHuman0412_NA__FLT3-5	tatatatcttgtggaaaggacgaaacACCGGCACGCTGTACAGACGACAAGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCGAGGCGCGCCGCTCCAGGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
2	HDAC1_chr1_32757816__NonTargetingControlGuideForHuman0412	HDAC1-1	gcgctcgcgcccggacgcgg	NonTargetingControlGuideForHuman0412_NA	GCACGCTGTACAGACGACAA	HDAC1_chr1_32757816	NonTargetingControlGuideForHuman0412	HDAC1-1__NonTargetingControlGuideForHuman0412_NA	tatatatcttgtggaaaggacgaaacACCGGACCGACTGACGGTAGGGACGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
10	FGFR2_chr10_123298215__NonTargetingControlGuideForHuman0412	FGFR2-13	gcgcggccgccACAAAGCTC	NonTargetingControlGuideForHuman0412_NA	GCACGCTGTACAGACGACAA	FGFR2_chr10_123298215	NonTargetingControlGuideForHuman0412	FGFR2-13__NonTargetingControlGuideForHuman0412_NA	tatatatcttgtggaaaggacgaaacACCGgcgcggccgccACAAAGCTCGTTTTgagacgTAGGGATAACAGGGTAATcgtctcGTTTGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT	146
""")
        # note that: list items are in sorted order by probe name, NonTargetingControlGuideForHuman0412 only occurs once
        # although it occurs multiple times in input, and the two seqs that weren't already all upper case have been
        # uppercased.
        expected_grna_names_and_seqs = [("FGFR2_chr10_123298215", "GCGCGGCCGCCACAAAGCTC"), # note: converted to upper
                                        ("FLT3_chr13_28636131", "GCGAGGCGCGCCGCTCCAGG"),
                                        ("HDAC1_chr1_32757816", "GCGCTCGCGCCCGGACGCGG"),  # note: converted to upper
                                        ("NonTargetingControlGuideForHuman0412", "GCACGCTGTACAGACGACAA")]

        # note that _read_in_construct_table is not being tested here; since this is a test of
        # _extract_grnas_from_construct_table, we're just assuming _read_in_construct_table works as expected
        #construct_table = _read_in_construct_table(construct_input, [1,2,6,3,4,7,5], rows_to_skip=1)

        #output_grna_names_and_seqs = _extract_unique_sets_across_a_and_b(construct_table)

        construct_names, output_grna_names_and_seqs = extract_construct_and_grna_info(construct_input, [1,2,6,3,4,7,5])

        self.assertEqual(len(expected_grna_names_and_seqs), len(output_grna_names_and_seqs))
        for x in range(0, len(expected_grna_names_and_seqs)):
            self.assertEqual(expected_grna_names_and_seqs[x], output_grna_names_and_seqs[x])

    # endregion

