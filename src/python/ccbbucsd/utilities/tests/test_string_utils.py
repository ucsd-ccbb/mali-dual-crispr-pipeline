# standard libraries
import unittest

# test library
import ccbbucsd.utilities.string_utils as ns_test


class TestFunctions(unittest.TestCase):
    def test_split_delimited_string_to_list_defaults(self):
        real_output = ns_test.split_delimited_string_to_list("1,3,5,9,blue")
        self.assertEqual(["1","3","5","9","blue"], real_output)

    def test_split_delimited_string_to_list_no_defaults(self):
        real_output = ns_test.split_delimited_string_to_list("1;3;5;9", int, ";")
        self.assertEqual([1,3,5,9], real_output)