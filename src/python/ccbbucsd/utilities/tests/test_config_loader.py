# standard libraries
import unittest

# test library
import ccbbucsd.utilities.config_loader as ns_test


class TestFunctions(unittest.TestCase):
    # No tests for load_config_parser_from_fp because its unique code just does file opening
    # No tests for get_default_config_fp because it is so simple

    # region load_config_parser
    def test_load_config_parser(self):
        input_config_str = """[DEFAULT]
main_dir: ${c4_2xlarge:main_dir}

[c4_2xlarge]
main_dir: /home/ec2-user/mali-dual-crispr-software
"""
        real_output = ns_test.load_config_parser(input_config_str)
        # the method under test returns a loaded config parser object, which has no useful tostring method.
        # I really don't want to check the type of the object.  I settle for checking that (a) the .get method works on
        # it, and (b) that one of the values that requires interpolation is correct, indicating the interpolation works
        real_interpolated_value = real_output.get(real_output.default_section, "main_dir")
        self.assertEqual("/home/ec2-user/mali-dual-crispr-software", real_interpolated_value)

    # end region

    def test_load_config_section_dict_default(self):
        input_config_str = """[DEFAULT]
main_dir: ${c4_2xlarge:main_dir}
data_dir: ${c4_2xlarge:data_dir}
num_processors: ${c4_2xlarge:num_processors}
keep_gzs: False
full_5p_r1: TATATATCTTGTGGAAAGGACGAAACACCG
len_of_seq_to_match = 19
num_allowed_mismatches = 1
code_dir: ${main_dir}/src/python
notebook_dir: ${main_dir}/notebooks
raw_data_dir: ${data_dir}/raw
interim_data_dir: ${data_dir}/interim
processed_data_dir: ${data_dir}/processed

[c4_2xlarge]
main_dir: /home/ec2-user/mali-dual-crispr-software
data_dir: /data
num_processors: 7
"""
        expected_output = {"main_dir":"/home/ec2-user/mali-dual-crispr-software",
                           "data_dir":"/data",
                           "keep_gzs":"False",
                           "full_5p_r1":"TATATATCTTGTGGAAAGGACGAAACACCG",
                           "len_of_seq_to_match":"19",
                           "num_allowed_mismatches":"1",
                           "code_dir":"/home/ec2-user/mali-dual-crispr-software/src/python",
                           "notebook_dir":"/home/ec2-user/mali-dual-crispr-software/notebooks",
                           "raw_data_dir":"/data/raw",
                           "interim_data_dir":"/data/interim",
                           "processed_data_dir":"/data/processed",
                           "num_processors":"7"
        }

        # Note: not testing load_config_parser here, just using it
        loaded_config_parser = ns_test.load_config_parser(input_config_str)

        real_output = ns_test.load_config_section_dict(loaded_config_parser)
        self.assertEqual(expected_output, real_output)


    def test_load_config_section_dict_not_default(self):
        input_config_str = """[DEFAULT]
main_dir: ${c4_2xlarge:main_dir}
data_dir: ${c4_2xlarge:data_dir}
num_processors: ${c4_2xlarge:num_processors}
keep_gzs: False
full_5p_r1: TATATATCTTGTGGAAAGGACGAAACACCG
len_of_seq_to_match = 19
num_allowed_mismatches = 1
code_dir: ${main_dir}/src/python
notebook_dir: ${main_dir}/notebooks
raw_data_dir: ${data_dir}/raw
interim_data_dir: ${data_dir}/interim
processed_data_dir: ${data_dir}/processed

[c4_2xlarge]
main_dir: /home/ec2-user/mali-dual-crispr-software
data_dir: /data
num_processors: 7

[laptop]
main_dir: /Users/Birmingham/Work/Repositories/ccbb_tickets_2017/mali-dual-crispr-software
data_dir: /Users/Birmingham/Work/Data
num_processors: 3
"""
        expected_output = {"main_dir":"/Users/Birmingham/Work/Repositories/ccbb_tickets_2017/mali-dual-crispr-software",
                           "data_dir":"/Users/Birmingham/Work/Data",
                           "keep_gzs": "False",
                           "full_5p_r1": "TATATATCTTGTGGAAAGGACGAAACACCG",
                           "len_of_seq_to_match": "19",
                           "num_allowed_mismatches": "1",
                           "code_dir": "/Users/Birmingham/Work/Repositories/ccbb_tickets_2017/mali-dual-crispr-software/src/python",
                           "notebook_dir": "/Users/Birmingham/Work/Repositories/ccbb_tickets_2017/mali-dual-crispr-software/notebooks",
                           "raw_data_dir": "/Users/Birmingham/Work/Data/raw",
                           "interim_data_dir": "/Users/Birmingham/Work/Data/interim",
                           "processed_data_dir": "/Users/Birmingham/Work/Data/processed",
                           "num_processors":"3"
        }

        # Note: not testing load_config_parser here, just using it
        loaded_config_parser = ns_test.load_config_parser(input_config_str)

        real_output = ns_test.load_config_section_dict(loaded_config_parser, "laptop")
        self.assertEqual(expected_output, real_output)