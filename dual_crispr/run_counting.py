# standard libraries
import argparse
import os

import ccbb_pyutils.config_loader as ns_config
import ccbb_pyutils.notebook_pipeliner as ns_pipeliner

from dual_crispr import dual_crispr_pipeliner as ns_dcpipe

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def _parse_cmd_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_name", help="short, alphanumeric human-readable name for the dataset to be analyzed")
    parser.add_argument("library_name", help="name of the construct library for the dataset to be analyzed")
    parser.add_argument("fastq_dir_path", help="absolute path to the folder in which the fastq data to be analyzed reside")
    parser.add_argument("output_dir_path", help="absolute path to folder in which new output directory should be created")
    parser.add_argument("--config", help="path to config file; if not specified, config file in default location will "
                                         "be used")
    args = parser.parse_args()
    return args.fastq_dir_path, args.output_dir_path, args.dataset_name, args.library_name, args.config


def _set_params(fastq_dir_path, outputs_dir_path, config_fp):
    len_of_seq_to_match_key = "len_of_seq_to_match"
    num_allowed_mismatches_key = "num_allowed_mismatches"

    # load the config file
    config_fp = ns_dcpipe.get_config_fp_or_default(config_fp)
    configparser = ns_config.load_config_parser_from_fp(config_fp)
    count_params = ns_config.load_config_section_dict(configparser, "count_pipeline")
    result = count_params.copy()

    # write the path params
    result[ns_dcpipe.DirectoryKeys.RAW_DATA.value] = fastq_dir_path
    result[ns_dcpipe.DirectoryKeys.INTERIM_DATA.value] = os.path.join(outputs_dir_path, "temporary_files")
    result[ns_dcpipe.DirectoryKeys.PROCESSED_DATA.value] = outputs_dir_path

    result["g_fastqs_dir"] = result[ns_dcpipe.DirectoryKeys.RAW_DATA.value]
    result["g_trimmed_fastqs_dir"] = result[ns_dcpipe.DirectoryKeys.INTERIM_DATA.value]
    result["g_filtered_fastqs_dir"] = result[ns_dcpipe.DirectoryKeys.INTERIM_DATA.value]

    # cast the params that aren't supposed to be strings
    result[len_of_seq_to_match_key] = int(result[len_of_seq_to_match_key])
    result[num_allowed_mismatches_key] = int(result[num_allowed_mismatches_key])

    return result


def main():
    fastq_dir_path, outputs_dir_path, dataset_name, library_name, config_fp = _parse_cmd_line_args()
    count_params = _set_params(fastq_dir_path, outputs_dir_path, config_fp)
    full_params = ns_dcpipe.generate_notebook_params(dataset_name, library_name, count_params, config_fp)
    # Note: second argument is the *function*, not results of calling the function
    ns_pipeliner.execute_run_from_full_params(full_params, ns_dcpipe.rename_param_names_as_global_vars)

if __name__ == '__main__':
    main()
