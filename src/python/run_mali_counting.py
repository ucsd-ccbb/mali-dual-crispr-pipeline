# standard libraries
import argparse
import os

# ccbb libraries
import ccbbucsd.utilities.analysis_run_prefixes as ns_runs
import ccbbucsd.utilities.config_loader as ns_config
import ccbbucsd.utilities.notebook_pipeliner as ns_pipeliner

# project-specific libraries
import dual_crispr_pipeliner as ns_dcpipe

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def _parse_cmd_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_dir_name", help="name of the folder in which the fastq data to be analyzed reside")
    parser.add_argument("dataset_name", help="short, alphanumeric human-readable name for the dataset to be analyzed")
    parser.add_argument("library_name", help="name of the construct library for the dataset to be analyzed")
    args = parser.parse_args()
    return args.fastq_dir_name, args.dataset_name, args.library_name


def _set_params(fastq_dir_name):
    # load the config file
    config_params = ns_dcpipe.get_machine_config_params()
    raw_dir = config_params[ns_dcpipe.DirectoryKeys.RAW_DATA.value]
    interim_dir = config_params[ns_dcpipe.DirectoryKeys.INTERIM_DATA.value]
    count_notebooks_list = config_params["count_notebooks"]

    result = {}
    result[ns_dcpipe.DirectoryKeys.RAW_DATA.value] = os.path.join(raw_dir, fastq_dir_name)
    result[ns_dcpipe.DirectoryKeys.INTERIM_DATA.value] = os.path.join(interim_dir, fastq_dir_name)
    result[ns_runs.get_notebook_names_list_key()] = count_notebooks_list

    return result


def main():
    fastq_dir_name, dataset_name, library_name = parse_cmd_line_args()
    count_params = set_params(fastq_dir_name)
    full_params = ns_dcpipe.generate_notebook_params(dataset_name, library_name, count_params)
    # Note: second argument is the *function*, not results of calling the function
    ns_pipeliner.execute_run_from_full_params(full_params, ns_dcpipe.rename_param_names_as_global_vars)

if __name__ == '__main__':
    main()
