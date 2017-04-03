# standard libraries
import io
import unittest

# project-specific libraries
from dual_crispr.count_files_and_dataframes import clip_count_header_suffix

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
    # region clip_count_header_suffix
    def test_clip_count_header_suffix_igm_clipped(self):
        header = "PGP1-MV4_d21_2_S9_trimmed53_len_filtered_counts"
        expected_output = "PGP1-MV4_d21_2"
        real_output = clip_count_header_suffix(header)
        self.assertEqual(expected_output, real_output)

    def test_clip_count_header_suffix_pipeline_clipped(self):
        header = "somekindaname_trimmed53_len_filtered_counts"
        expected_output = "somekindaname"
        real_output = clip_count_header_suffix(header)
        self.assertEqual(expected_output, real_output)

    def test_clip_count_header_suffix_unclipped(self):
        header = "PGP1-MV4_d21_2_S946"
        real_output = clip_count_header_suffix(header)
        self.assertEqual(header, real_output)

    # end region
