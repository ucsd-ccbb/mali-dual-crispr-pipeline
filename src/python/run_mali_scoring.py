# standard libraries
import argparse
import warnings

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
    # examples:
    # human_readable_name = 20160627HeLaA549
    # library_name = CV4
    # day_timepoints_str = 3,14,20,28
    # --test

    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_name", help="short, alphanumeric human-readable name for the dataset to be analyzed")
    parser.add_argument("library_name", help="name of the construct library for the dataset to be analyzed")
    parser.add_argument("count_fps_or_dirs", help="a comma-separated list of the paths to the file of counts to be "
                                                  "scored, or to the directories in which the *_combined_counts.txt "
                                                  "file to be scored resides")
    parser.add_argument("day_timepoints_str", help="a comma-separated list containing the time points (in order) at "
                                                   "which data were collected")
    parser.add_argument("--test", help="run with set seed and only two iterations, suitable for testing ONLY",
                        action="store_true")
    parser.add_argument("--config", help="path to config file; if not specified, config file in default location will "
                                         "be used")
    args = parser.parse_args()
    return args.dataset_name, args.library_name, args.counts_fps_or_dirs, args.day_timepoints_str, args.test, args.config


def _set_params(count_fps_or_dirs, day_timepoints_str, is_test, config_fp):
    test_config_section_key = "test"
    score_notebooks_key = "score_notebooks"
    time_prefixes_key = "time_prefixes"
    count_fps_or_dirs_key = "count_fps_or_dirs"
    day_timepoints_str_key = "day_timepoints_str"
    use_seed_key = "use_seed"
    num_iterations_key = "num_iterations"

    # load the config file
    config_parser = ns_config.load_config_parser_from_fp(config_fp)
    # delimited string NOT split here--done in shared code in dual_crispr_pipeliner
    score_notebooks_list_str = config_parser.get(config_parser.default_section, score_notebooks_key)

    # Note: the time_prefixes_str comma-delimited string value below need to be converted to a list, but the conversion
    # is NOT being done here--it is done in the notebook, because if users run the notebook directly,
    # they will have to put in a comma-delimited string there, so the notebook needs to know how to deal with it.
    time_prefixes_str = config_parser.get(config_parser.default_section, time_prefixes_key)

    # Inputs below *are* converted here because notebook users can directly input booleans and ints into
    # the notebook params (unlike lists, which nbparameterise won't accept)
    score_config_section = test_config_section_key if is_test else config_parser.default_section
    use_seed = config_parser.getboolean(score_config_section, use_seed_key)
    num_iterations = config_parser.getint(score_config_section, num_iterations_key)

    # if test parameters are detected, remind user! Results should not be used for real analysis
    if use_seed or score_config_section == test_config_section_key:
        warnings.warn('Scoring is running in TEST MODE; do not use results for data analysis!')

    result = {}
    result[ns_runs.get_notebook_names_list_key()] = score_notebooks_list_str
    result[time_prefixes_key] = time_prefixes_str
    result[count_fps_or_dirs_key] = count_fps_or_dirs
    result[day_timepoints_str_key] = day_timepoints_str
    result[use_seed_key] = use_seed
    result[num_iterations_key] = num_iterations

    return result


def main():
    dataset_name, library_name, counts_fps_or_dirs, day_timepoints_str, is_test, config_fp = _parse_cmd_line_args()
    score_params = _set_params(counts_fps_or_dirs, day_timepoints_str, is_test, config_fp)
    full_params = ns_dcpipe.generate_notebook_params(dataset_name, library_name, score_params)
    # Note: second argument is the *function*, not results of calling the function
    ns_pipeliner.execute_run_from_full_params(full_params, ns_dcpipe.rename_param_names_as_global_vars)


if __name__ == '__main__':
    main()