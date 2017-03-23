# standard libraries
import os
import tempfile
import unittest

import ccbbucsd.utilities.config_loader as ns_config

# test library
import dual_crispr_pipeliner as ns_test


class TestFunctions(unittest.TestCase):
    # region _validate_expt_name
    def test__validate_expt_name_invalid(self):
        with self.assertRaises(ValueError):
            ns_test._validate_expt_name("PGP1_MV4")

    def test__validate_expt_name_valid(self):
        ns_test._validate_expt_name("PGP1MV4")
        self.assertTrue(True)  # if it got this far, the test passed

    # end region

    def test_rename_param_names_as_global_vars(self):
        input_dict = {"x":5,
                      "g_y":"blue"}
        expected_output = {"g_x":5,
                           "g_y":"blue"}
        real_output = ns_test.rename_param_names_as_global_vars(input_dict)
        self.assertEqual(expected_output, real_output)

    def test__verify_or_make_data_subdirs(self):
        tempdir = tempfile.TemporaryDirectory()
        input_dict = {"raw_data_dir": os.path.join(tempdir.name, "raw_test"),
                      "interim_data_dir": os.path.join(tempdir.name, "in_test"),
                      "processed_data_dir": os.path.join(tempdir.name, "pro_test")}
        for curr_key, curr_val in input_dict.items():
            self.assertFalse(os.path.exists(curr_val))

        ns_test._verify_or_make_data_subdirs(input_dict)

        for curr_key, curr_val in input_dict.items():
            self.assertTrue(os.path.exists(curr_val))

    def test__get_machine_config_params(self):
        config_string = """[DEFAULT]
machine_configuration = c4_2xlarge
keep_gzs: False

[laptop]
main_dir: /Users/Birmingham/Work/Repositories/ccbb_tickets_2017/
data_dir: /Users/Birmingham/Work/Data
num_processors: 3
code_dir: ${main_dir}/src/python

[c4_2xlarge]
main_dir: /home/ec2-user
data_dir: /data
num_processors: 7
code_dir: ${main_dir}/src/python
"""

        expected_output = {"machine_configuration": "c4_2xlarge",
                           "keep_gzs": "False",
                           "main_dir": "/home/ec2-user",
                           "data_dir": "/data",
                           "num_processors": "7",
                           "code_dir": "/home/ec2-user/src/python"
                           }

        configparser = ns_config.load_config_parser(config_string)
        real_output = ns_test._get_machine_config_params(configparser)
        self.assertEqual(expected_output, real_output)

    def test_set_up_data_subdirs_and_get_machine_configs(self):
        tempdir = tempfile.TemporaryDirectory()

        config_string = """[DEFAULT]
machine_configuration = laptop
keep_gzs: False

[laptop]
main_dir: REPLACE
data_dir: REPLACE
num_processors: 3
raw_data_dir: ${data_dir}/raw
interim_data_dir: ${data_dir}/interim
processed_data_dir: ${data_dir}/processed

[c4_2xlarge]
main_dir: /home/ec2-user
data_dir: /data
num_processors: 7
raw_data_dir: ${data_dir}/raw
interim_data_dir: ${data_dir}/interim
processed_data_dir: ${data_dir}/processed
"""

        config_string = config_string.replace("REPLACE", tempdir.name)
        configparser = ns_config.load_config_parser(config_string)

        expected_output = {"machine_configuration": "laptop",
                           "keep_gzs": "False",
                           "main_dir": tempdir.name,
                           "data_dir": tempdir.name,
                           "num_processors": "3"}
        dirs_dict = {ns_test.DirectoryKeys.RAW_DATA.value: os.path.join(tempdir.name, "raw"),
                      ns_test.DirectoryKeys.INTERIM_DATA.value: os.path.join(tempdir.name, "interim"),
                      ns_test.DirectoryKeys.PROCESSED_DATA.value: os.path.join(tempdir.name, "processed")}
        for curr_key in dirs_dict:
            self.assertFalse(os.path.exists(dirs_dict[curr_key]))

        real_output = ns_test.set_up_data_subdirs_and_get_machine_configs(configparser)

        for curr_key, curr_val in dirs_dict.items():
            self.assertTrue(os.path.exists(curr_val))

        expected_output.update(dirs_dict)
        self.assertEqual(expected_output, real_output)

    def test_generate_notebook_params(self):
        tempdir = tempfile.TemporaryDirectory()

        config_string = """[DEFAULT]
machine_configuration = laptop
keep_gzs: False

[laptop]
main_dir: REPLACE
data_dir: REPLACE
num_processors: 3
code_dir: ${main_dir}/src/python
notebook_dir: ${main_dir}/notebooks
libraries_dir: ${main_dir}/library_definitions
raw_data_dir: ${data_dir}/raw
interim_data_dir: ${data_dir}/interim
processed_data_dir: ${data_dir}/processed

[c4_2xlarge]
main_dir: /home/ec2-user
data_dir: /data
num_processors: 7
code_dir: ${main_dir}/src/python
notebook_dir: ${main_dir}/notebooks
libraries_dir: ${main_dir}/library_definitions
raw_data_dir: ${data_dir}/raw
interim_data_dir: ${data_dir}/interim
processed_data_dir: ${data_dir}/processed
"""

        library_file_contents = """# library_name = CV4
# min_trimmed_grna_len = 19
# max_trimmed_grna_len = 21
construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq	FinalSequence	Gene_A_chr	Gene_A_pos	Gene_B_chr	Gene_B_pos
BRCA1_chr17_41276018__NonTargetingControlGuideForHuman0412	BRCA1	BRCA1_chr17_41276018	TCTTGTGCTGACTTACCAGA	NonTargetingControlGuideForHuman0412	NonTargetingControlGuideForHuman0412	GCACGCTGTACAGACGACAA	tatatatcttgtggaaaggacgaaacACCGTCTTGTGCTGACTTACCAGAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGGCACGCTGTACAGACGACAAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	chr17	41276018	NA	NA
NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972	NonTargetingControlGuideForHuman0352	NonTargetingControlGuideForHuman0352	GCCATTCTAGTCCCGGCATA	SETD2	SETD2_chr3_47142972	AACGTTAACTCTGAGCCTGA	tatatatcttgtggaaaggacgaaacACCGGCCATTCTAGTCCCGGCATAGTTTcAGAGCTAtgctgGAAActgcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAGtCGCCCaGTCTCAGATAGATCCGACGCCGCCATCTCTAGGCCCGCGCCGGCCCCCTCGCACAGACTTGTGGGAGAAGCTCGGCTACTCCCCTGCCCCGGTTAATTTGCATATAATATTTCCTAGTAACTATAGAGGCTTAATGTGCGATAAAAGAGAACGTTAACTCTGAGCCTGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG	NA	NA	chr3	47142972
"""

        first_nb = """{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "A non-code cell"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "x = 5\\n",
    "g_prepped_counts_run_prefix = \\"blue\\"\\n",
    "g_prepped_counts_dir = \\"red\\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": false
   },
   "outputs": [],
   "source": [
    "print(\\"x is {0} and g_prepped_counts_run_prefix is {1}\\".format(x, g_prepped_counts_run_prefix))"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.4.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 1
}
"""

        second_nb = """{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "A non-code cell"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "x = 5\\n",
    "g_micro_run_prefix = \\"blue\\"\\n",
    "g_macro_dir = \\"red\\"\\n",
    "g_macro_run_prefix = \\"puce\\"\\n",
    "y = \\"cerulian\\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": false
   },
   "outputs": [],
   "source": [
    "print(\\"x is {0} and y is {1}\\".format(x, y))"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.4.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 1
}
"""
        # set up the mock config parser
        config_string = config_string.replace("REPLACE", tempdir.name)
        configparser = ns_config.load_config_parser(config_string)

        # set up the mock library file
        temp_library_dir = configparser.get("laptop", "libraries_dir")
        os.mkdir(temp_library_dir)
        temp_library_fp = os.path.join(configparser.get("laptop", "libraries_dir"), "temp_lib.txt")
        with open(temp_library_fp, "w") as f:
            f.write(library_file_contents)

        # set up the mock notebooks
        temp_nb_dir = configparser.get("laptop", "notebook_dir")
        os.mkdir(temp_nb_dir)
        temp_nb_fp = os.path.join(temp_nb_dir, "first_nb.ipynb")
        with open(temp_nb_fp, "w") as f:
            f.write(first_nb)
        temp_nb_fp = os.path.join(temp_nb_dir, "second_nb.ipynb")
        with open(temp_nb_fp, "w") as f:
            f.write(second_nb)

        # set up the mock argument-based parameters
        arg_based_params_input = {"interim_data_dir":"an overwrite",
                                  "notebook_basenames_list": "first_nb.ipynb, second_nb.ipynb"}

        real_output = ns_test.generate_notebook_params("test", "CV4", arg_based_params_input, configparser)

        # check that default values from the config "file" were added
        self.assertEqual("laptop", real_output["machine_configuration"])
        self.assertEqual("False", real_output["keep_gzs"])

        # check that params from the selected machine config in the config "file" were added
        self.assertEqual("3", real_output["num_processors"])
        # testing here that all dir keys exist, but not what they are since they
        # are being generated on the fly in the temp directory
        for curr_key in ['main_dir','data_dir', "raw_data_dir", 'processed_data_dir',
                         'code_dir', 'libraries_dir', 'notebook_dir']:  # note interim_dir tested later
            self.assertTrue(curr_key in real_output)

        # check that the parameters from the library file were added
        self.assertEqual("CV4", real_output["library_name"])
        self.assertEqual("19", real_output["min_trimmed_grna_len"])
        self.assertEqual("21", real_output["max_trimmed_grna_len"])
        self.assertEqual(temp_library_fp, real_output["library_fp"])

        # check that the run prefix and run dir were added
        self.assertTrue("run_prefix" in real_output)
        self.assertTrue("run_dir" in real_output)
        run_prefix = real_output["run_prefix"]
        run_dir = real_output["run_dir"]

        # check that the comma-delimited notebook names string that was put in
        # has been converted to an actual list
        self.assertEqual(['first_nb.ipynb', 'second_nb.ipynb'], real_output['notebook_basenames_list'])

        # check that run prefix and dir variables were found in notebooks and set to
        # the values for this run
        self.assertEqual(run_prefix, real_output["g_prepped_counts_run_prefix"])
        self.assertEqual(run_prefix, real_output["g_micro_run_prefix"])
        self.assertEqual(run_prefix, real_output["g_macro_run_prefix"])
        self.assertEqual(run_dir, real_output["g_prepped_counts_dir"])
        self.assertEqual(run_dir, real_output["g_macro_dir"])

        # check that when an arg-based param conflicts with a config-based param,
        # the arg-based param overwrites the config-based one
        self.assertEqual("an overwrite", real_output["interim_data_dir"])

