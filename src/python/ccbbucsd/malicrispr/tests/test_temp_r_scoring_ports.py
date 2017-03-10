# standard libraries
import io
import unittest

# third-party libraries
import pandas

# project-specific libraries
from ccbbucsd.malicrispr.temp_r_scoring_ports import get_subset_of_df_rows

ABUNDANCE_THRESHS = [-19.0, -18.5, -18.5, -19.0, -19.0, -19.0, -19.0, -19.0]
REP1_ABUND_THRESHS = [-19.0, -18.5, -19.0, -19.0]
REP2_ABUND_THRESHS = [-18.5, -19.0, -19.0, -19.0]


class TestFunctions(unittest.TestCase):
    # region _clip_count_header_suffix
    def test_get_subset_of_df_rows_include_by_colnames(self):

        input_df_str = io.StringIO("""col1;col2;col3
1;4.4;99
2;4.5;200
3;4.7;65
4;3.2;140
""")

        expected_output_df_str = io.StringIO("""col1;col2;col3
2;4.5;200
4;3.2;140
""")

        input_df = pandas.read_csv(input_df_str, sep=";")
        expected_outout_df = pandas.read_csv(expected_output_df_str, sep=";")

        output_df = get_subset_of_df_rows(input_df, [140, 200], "col3", False)
        self.assertEqual(expected_outout_df, output_df)