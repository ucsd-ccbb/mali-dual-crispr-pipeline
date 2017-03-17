# standard libraries
import io
import unittest
import warnings

# project-specific libraries
from ccbbucsd.malicrispr.constructs_settings_extracter import _validate_settings_values, \
    _validate_settings_keys, _LIBRARY_NAME_KEY, _LIBRARY_FP_KEY, _LIBRARY_SETTINGS_KEYS, \
    _validate_library_id, _id_and_validate_library, _id_and_trim_settings_line, _mine_library_file

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
    # id_library_info not tested because it just calls _extract_library_info_from_files and _id_and_validate_library
    # _extract_library_info_from_files not tested because its unique code is file-system and file-opening

    # region _validate_settings_values
    def test__validate_settings_values_valid(self):
        input_dict = {"key1":"blusky",
                      "key2":"goldish",
                      "key3":"reddish"}

        _validate_settings_values("CV4", input_dict)

    def test__validate_settings_values_error(self):
        input_dict = {"key1":"blusky",
                      "key2":"",
                      "key3":"reddish"}

        with self.assertRaises(ValueError):
            _validate_settings_values("CV4", input_dict)

    # endregion

    # region _validate_settings_keys
    def test__validate_settings_keys_valid(self):
        input_settings = {_LIBRARY_NAME_KEY: "CV4",
                             _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                             _LIBRARY_SETTINGS_KEYS[0]: "19",
                             _LIBRARY_SETTINGS_KEYS[1]: "21"}

        _validate_settings_keys("CV4", input_settings)

    def test__validate_settings_keys_valid_w_warning(self):
        input_settings = {_LIBRARY_NAME_KEY: "CV4",
                         _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                         _LIBRARY_SETTINGS_KEYS[0]: "19",
                         _LIBRARY_SETTINGS_KEYS[1]: "21",
                         "key2":"goldish"}

        with warnings.catch_warnings(record=True) as warnings_list:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")

            _validate_settings_keys("CV4", input_settings)

            assert len(warnings_list) == 1
            assert "ignored" in str(warnings_list[-1].message)

    def test__validate_settings_keys_error(self):
        expected_settings = {_LIBRARY_NAME_KEY: "CV4",
                             _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                             _LIBRARY_SETTINGS_KEYS[0]: "19",
                             _LIBRARY_SETTINGS_KEYS[1]: "21"}
        input_settings = {_LIBRARY_NAME_KEY: "CV4",
                        _LIBRARY_FP_KEY: "/some/path/to/CV4.txt"}

        with self.assertRaises(ValueError):
            _validate_settings_keys("CV4", input_settings)

    def test__validate_settings_keys_error_w_warning(self):
        # inputs are missing the expected key "library_name" but have a typo of it instead
        input_settings = {"libraryname": "CV4",
                          _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                          _LIBRARY_SETTINGS_KEYS[0]: "19",
                          _LIBRARY_SETTINGS_KEYS[1]: "21"}

        with warnings.catch_warnings(record=True) as warnings_list:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")

            with self.assertRaises(ValueError):
                _validate_settings_keys("CV4", input_settings)

            assert len(warnings_list) == 1
            assert "ignored" in str(warnings_list[-1].message)

    # endregion

    # region _validate_library_id
    def test__validate_library_id_valid(self):
        input_dicts_list = [
            {_LIBRARY_NAME_KEY: "CV4",
             _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
             _LIBRARY_SETTINGS_KEYS[0]: "19",
             _LIBRARY_SETTINGS_KEYS[1]: "21"}
        ]

        _validate_library_id("CV4", input_dicts_list)

    def test__validate_library_id_error_too_few(self):
        with self.assertRaises(ValueError):
            _validate_library_id("CV4", [])

    def test__validate_library_id_error_too_many(self):
        input_dicts_list = [
            {_LIBRARY_NAME_KEY: "CV4",
             _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
             _LIBRARY_SETTINGS_KEYS[0]: "19",
             _LIBRARY_SETTINGS_KEYS[1]: "21"},
            {_LIBRARY_NAME_KEY: "CV4",
            _LIBRARY_FP_KEY: "/old/path/to/CV4.txt",
            _LIBRARY_SETTINGS_KEYS[0]: "19",
            _LIBRARY_SETTINGS_KEYS[1]: "21"}
        ]

        with self.assertRaises(ValueError):
            _validate_library_id("CV4", input_dicts_list)

    # endregion

    # region _id_and_validate_library
    def test__id_and_validate_library_valid(self):
        expected_output = {_LIBRARY_NAME_KEY: "CV4",
             _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
             _LIBRARY_SETTINGS_KEYS[0]: "19",
             _LIBRARY_SETTINGS_KEYS[1]: "21"}

        input_dicts_list = [expected_output]

        real_output = _id_and_validate_library("CV4", input_dicts_list)
        self.assertEqual(expected_output, real_output)

    def test__id_and_validate_library_too_many_libs(self):
        expected_output = {_LIBRARY_NAME_KEY: "CV4",
                           _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                           _LIBRARY_SETTINGS_KEYS[0]: "19",
                           _LIBRARY_SETTINGS_KEYS[1]: "21"}

        input_dicts_list = [expected_output,
                            {_LIBRARY_NAME_KEY: "CV4",
                             _LIBRARY_FP_KEY: "/old/path/to/CV4.txt",
                             _LIBRARY_SETTINGS_KEYS[0]: "19",
                             _LIBRARY_SETTINGS_KEYS[1]: "21"}
                            ]

        with self.assertRaises(ValueError):
            _id_and_validate_library("CV4", input_dicts_list)

    def test__id_and_validate_library_too_few_libs(self):
        with self.assertRaises(ValueError):
            _id_and_validate_library("CV4", [])

    # Note: Not testing all the different combinations of raising of warnings for ignored keys here because I already
    # did so above and it isn't core functionality, so it doesn't seem worth the extra effort to write more tests for it

    def test__id_and_validate_library_missing_and_ignored_keys(self):
        # inputs are missing the expected key "library_name" but have a typo of it instead
        input_dicts_list = [{"libraryname": "CV4",
                          _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                          _LIBRARY_SETTINGS_KEYS[0]: "19",
                          _LIBRARY_SETTINGS_KEYS[1]: "21"}]

        with warnings.catch_warnings(record=True) as warnings_list:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")

            with self.assertRaises(ValueError):
                _id_and_validate_library("CV4", input_dicts_list)

            assert len(warnings_list) == 1
            assert "ignored" in str(warnings_list[-1].message)

    def test__id_and_validate_library_missing_values(self):
        input_dicts_list = [{_LIBRARY_NAME_KEY: "CV4",
                             _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                             _LIBRARY_SETTINGS_KEYS[0]: "",
                             _LIBRARY_SETTINGS_KEYS[1]: "21"}]

        with self.assertRaises(ValueError):
            _id_and_validate_library("CV4", input_dicts_list)

    # endregion

    # region _id_and_trim_settings_line
    def test__id_and_trim_settings_line_valid(self):
        expected_output = ('library_name', 'CV4')
        real_output = _id_and_trim_settings_line("#library_name=CV4")
        self.assertEqual(expected_output, real_output)

        real_output2 = _id_and_trim_settings_line("# library_name = CV4       ")
        self.assertEqual(expected_output, real_output2)

    def test__id_and_trim_settings_line_none_no_comment_char(self):
        real_output = _id_and_trim_settings_line("library_name=CV4")
        self.assertIsNone(real_output)

    def test__id_and_trim_settings_line_none_no_split(self):
        real_output = _id_and_trim_settings_line("# library_name: CV4")
        self.assertIsNone(real_output)

    def test__id_and_trim_settings_line_none_multi_split(self):
        real_output = _id_and_trim_settings_line("# library_name=CV4=old")
        self.assertIsNone(real_output)

    def test__id_and_trim_settings_line_none_one_sided_split(self):
        real_output = _id_and_trim_settings_line("# library_name = ")
        self.assertIsNone(real_output)

    # endregion

    # region _mine_library_file
    def test__mine_library_file_valid(self):
        file_contents_str = """# library_name = CV4
# min_trimmed_grna_len = 19
# max_trimmed_grna_len = 21
construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq	FinalSequence	Gene_A_chr	Gene_A_pos	Gene_B_chr	Gene_B_pos
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	TCTTGTGCTGACTTACCAGA	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	tatatatcttgtggaaaggacgaaacACCGTCTTGTGCTGACTTACCAGAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	chr17	41276018	NA	NA
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	NonTargetingControlGuideForHuman0352	GCCATTCTAGTCCCGGCATA	SETD2	SETD2_chr3_47142972	AACGTTAACTCTGAGCCTGA	tatatatcttgtggaaaggacgaaacACCGGCCATTCTAGTCCCGGCATAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGAACGTTAACTCTGAGCCTGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	NA	NA	chr3	47142972
"""
        expected_output = {_LIBRARY_NAME_KEY: "CV4",
                           _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                           _LIBRARY_SETTINGS_KEYS[1]: "19",
                           _LIBRARY_SETTINGS_KEYS[0]: "21"}

        file_contents_obj = io.StringIO(file_contents_str)
        real_output = _mine_library_file("CV4", "/some/path/to/CV4.txt", file_contents_obj)
        self.assertEqual(expected_output, real_output)

    def test__mine_library_file_none_no_comment(self):
        file_contents_str = """library_name = CV4
# min_trimmed_grna_len = 19
# max_trimmed_grna_len = 21
construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq	FinalSequence	Gene_A_chr	Gene_A_pos	Gene_B_chr	Gene_B_pos
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	TCTTGTGCTGACTTACCAGA	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	tatatatcttgtggaaaggacgaaacACCGTCTTGTGCTGACTTACCAGAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	chr17	41276018	NA	NA
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	NonTargetingControlGuideForHuman0352	GCCATTCTAGTCCCGGCATA	SETD2	SETD2_chr3_47142972	AACGTTAACTCTGAGCCTGA	tatatatcttgtggaaaggacgaaacACCGGCCATTCTAGTCCCGGCATAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGAACGTTAACTCTGAGCCTGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	NA	NA	chr3	47142972
"""

        file_contents_obj = io.StringIO(file_contents_str)
        real_output = _mine_library_file("CV4", "/some/path/to/CV4.txt", file_contents_obj)
        # first line doesn't start with a #, so file is not recognized as a potential library file
        self.assertIsNone(real_output)

    def test__mine_library_file_none_first_header_wrong(self):
        file_contents_str = """# libraryname = CV4
# min_trimmed_grna_len = 19
# max_trimmed_grna_len = 21
construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq	FinalSequence	Gene_A_chr	Gene_A_pos	Gene_B_chr	Gene_B_pos
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	TCTTGTGCTGACTTACCAGA	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	tatatatcttgtggaaaggacgaaacACCGTCTTGTGCTGACTTACCAGAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	chr17	41276018	NA	NA
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	NonTargetingControlGuideForHuman0352	GCCATTCTAGTCCCGGCATA	SETD2	SETD2_chr3_47142972	AACGTTAACTCTGAGCCTGA	tatatatcttgtggaaaggacgaaacACCGGCCATTCTAGTCCCGGCATAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGAACGTTAACTCTGAGCCTGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	NA	NA	chr3	47142972
"""

        file_contents_obj = io.StringIO(file_contents_str)
        real_output = _mine_library_file("CV4", "/some/path/to/CV4.txt", file_contents_obj)
        # first line doesn't have expected key (library_name) after comment, so file is not recognized as a potential
        # library file
        self.assertIsNone(real_output)

    def test__mine_library_file_none_not_desired_library(self):
        file_contents_str = """# library_name = MV4
# min_trimmed_grna_len = 19
# max_trimmed_grna_len = 21
construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq	FinalSequence	Gene_A_chr	Gene_A_pos	Gene_B_chr	Gene_B_pos
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	TCTTGTGCTGACTTACCAGA	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	tatatatcttgtggaaaggacgaaacACCGTCTTGTGCTGACTTACCAGAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	chr17	41276018	NA	NA
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	NonTargetingControlGuideForHuman0352	GCCATTCTAGTCCCGGCATA	SETD2	SETD2_chr3_47142972	AACGTTAACTCTGAGCCTGA	tatatatcttgtggaaaggacgaaacACCGGCCATTCTAGTCCCGGCATAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGAACGTTAACTCTGAGCCTGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	NA	NA	chr3	47142972
"""

        file_contents_obj = io.StringIO(file_contents_str)
        real_output = _mine_library_file("CV4", "/some/path/to/CV4.txt", file_contents_obj)
        # first line indicates this isn't the library we want, so nothing returned
        self.assertIsNone(real_output)

    def test__mine_library_file_partial_contents_missing_key(self):
        file_contents_str = """# library_name = CV4
# min_trimmed_grna_len = 19
construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq	FinalSequence	Gene_A_chr	Gene_A_pos	Gene_B_chr	Gene_B_pos
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	TCTTGTGCTGACTTACCAGA	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	tatatatcttgtggaaaggacgaaacACCGTCTTGTGCTGACTTACCAGAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	chr17	41276018	NA	NA
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	NonTargetingControlGuideForHuman0352	GCCATTCTAGTCCCGGCATA	SETD2	SETD2_chr3_47142972	AACGTTAACTCTGAGCCTGA	tatatatcttgtggaaaggacgaaacACCGGCCATTCTAGTCCCGGCATAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGAACGTTAACTCTGAGCCTGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	NA	NA	chr3	47142972
"""

        expected_output = {_LIBRARY_NAME_KEY: "CV4",
                           _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                           _LIBRARY_SETTINGS_KEYS[1]: "19"}

        file_contents_obj = io.StringIO(file_contents_str)
        real_output = _mine_library_file("CV4", "/some/path/to/CV4.txt", file_contents_obj)
        # one of required keys is missing from file, so output dict will be one entry short
        self.assertEqual(expected_output, real_output)

    def test__mine_library_file_partial_contents_malformed_key(self):
            file_contents_str = """# library_name = CV4
# min_trimmed_grna_len = 19
# max_trimmed_grna_len: 21
construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq	FinalSequence	Gene_A_chr	Gene_A_pos	Gene_B_chr	Gene_B_pos
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	TCTTGTGCTGACTTACCAGA	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	tatatatcttgtggaaaggacgaaacACCGTCTTGTGCTGACTTACCAGAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	chr17	41276018	NA	NA
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	NonTargetingControlGuideForHuman0352	GCCATTCTAGTCCCGGCATA	SETD2	SETD2_chr3_47142972	AACGTTAACTCTGAGCCTGA	tatatatcttgtggaaaggacgaaacACCGGCCATTCTAGTCCCGGCATAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGAACGTTAACTCTGAGCCTGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	NA	NA	chr3	47142972
"""

            expected_output = {_LIBRARY_NAME_KEY: "CV4",
                               _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                               _LIBRARY_SETTINGS_KEYS[1]: "19"}

            file_contents_obj = io.StringIO(file_contents_str)
            real_output = _mine_library_file("CV4", "/some/path/to/CV4.txt", file_contents_obj)
            # one of required keys is malformed in the file, so output dict will be one entry short
            self.assertEqual(expected_output, real_output)

        # end region
