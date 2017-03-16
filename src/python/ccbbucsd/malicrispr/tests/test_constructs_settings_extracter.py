# standard libraries
import io
import unittest
import warnings

# project-specific libraries
from ccbbucsd.malicrispr.constructs_settings_extracter import _validate_settings_values, \
    _validate_settings_keys, _LIBRARY_NAME_KEY, _LIBRARY_FP_KEY, _LIBRARY_SETTINGS_KEYS, \
    _validate_library_id

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class TestFunctions(unittest.TestCase):
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
                             _LIBRARY_SETTINGS_KEYS[0]: 19,
                             _LIBRARY_SETTINGS_KEYS[1]: 21}

        _validate_settings_keys("CV4", input_settings)

    def test__validate_settings_keys_valid_w_warning(self):
        input_settings = {_LIBRARY_NAME_KEY: "CV4",
                         _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
                         _LIBRARY_SETTINGS_KEYS[0]: 19,
                         _LIBRARY_SETTINGS_KEYS[1]: 21,
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
                             _LIBRARY_SETTINGS_KEYS[0]: 19,
                             _LIBRARY_SETTINGS_KEYS[1]: 21}
        input_settings = {_LIBRARY_NAME_KEY: "CV4",
                        _LIBRARY_FP_KEY: "/some/path/to/CV4.txt"}

        with self.assertRaises(ValueError):
            _validate_settings_keys("CV4", input_settings)

    # endregion

    # region _validate_library_id
    def test__validate_library_id_valid(self):
        input_dicts_list = [
            {_LIBRARY_NAME_KEY: "CV4",
             _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
             _LIBRARY_SETTINGS_KEYS[0]: 19,
             _LIBRARY_SETTINGS_KEYS[1]: 21}
        ]

        _validate_library_id("CV4", input_dicts_list)

    def test__validate_library_id_error_too_few(self):
        with self.assertRaises(ValueError):
            _validate_library_id("CV4", [])

    def test__validate_library_id_error_too_many(self):
        input_dicts_list = [
            {_LIBRARY_NAME_KEY: "CV4",
             _LIBRARY_FP_KEY: "/some/path/to/CV4.txt",
             _LIBRARY_SETTINGS_KEYS[0]: 19,
             _LIBRARY_SETTINGS_KEYS[1]: 21},
            {_LIBRARY_NAME_KEY: "CV4",
            _LIBRARY_FP_KEY: "/old/path/to/CV4.txt",
            _LIBRARY_SETTINGS_KEYS[0]: 19,
            _LIBRARY_SETTINGS_KEYS[1]: 21}
        ]

        with self.assertRaises(ValueError):
            _validate_library_id("CV4", input_dicts_list)

    # endregion
