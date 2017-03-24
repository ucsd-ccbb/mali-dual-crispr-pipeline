# standard libraries
import argparse
import distutils.util
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
    parser.add_argument("--config", help="path to config file; if not specified, config file in default location will "
                                         "be used")
    args = parser.parse_args()
    return args.fastq_dir_name, args.dataset_name, args.library_name, args.config


def _set_params(fastq_dir_name, config_fp):
    keep_gz_key = "keep_gzs"
    len_of_seq_to_match_key = "len_of_seq_to_match"
    num_allowed_mismatches_key = "num_allowed_mismatches"


    # load the config file
    config_params = ns_dcpipe.get_machine_config_params(config_fp)
    raw_dir = config_params[ns_dcpipe.DirectoryKeys.RAW_DATA.value]
    interim_dir = config_params[ns_dcpipe.DirectoryKeys.INTERIM_DATA.value]

    configparser = ns_config.load_config_parser_from_fp(config_fp)
    count_params = ns_config.load_config_section_dict(configparser, "count_pipeline")

    result = count_params.copy()
    # overwrite the params that need to be more specific
    result[ns_dcpipe.DirectoryKeys.RAW_DATA.value] = os.path.join(raw_dir, fastq_dir_name)
    result[ns_dcpipe.DirectoryKeys.INTERIM_DATA.value] = os.path.join(interim_dir, fastq_dir_name)

    # cast the params that aren't supposed to be strings
    result[keep_gz_key] = bool(distutils.util.strtobool(result[keep_gz_key]))
    result[len_of_seq_to_match_key] = int(result[len_of_seq_to_match_key])
    result[num_allowed_mismatches_key] = int(result[num_allowed_mismatches_key])

    return result


def main():
    fastq_dir_name, dataset_name, library_name, config_fp = _parse_cmd_line_args()
    count_params = _set_params(fastq_dir_name, config_fp)
    full_params = ns_dcpipe.generate_notebook_params(dataset_name, library_name, count_params, config_fp)
    # Note: second argument is the *function*, not results of calling the function
    ns_pipeliner.execute_run_from_full_params(full_params, ns_dcpipe.rename_param_names_as_global_vars)

if __name__ == '__main__':
    main()
