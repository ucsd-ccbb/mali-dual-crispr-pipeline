# standard libraries
import argparse
import distutils.util
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
    return args.dataset_name, args.library_name, args.count_fps_or_dirs, args.day_timepoints_str, args.test, args.config


def _set_params(count_fps_or_dirs, day_timepoints_str, is_test, config_fp):
    test_config_section_key = "test"
    count_fps_or_dirs_key = "count_fps_or_dirs"
    min_count_limit_key = "min_count_limit"
    max_fraction_acceptable_spline_density_diff_key = "max_fraction_acceptable_spline_density_diff"
    max_fraction_counts_excluded_key = "max_fraction_counts_excluded"
    day_timepoints_str_key = "day_timepoints_str"
    use_seed_key = "use_seed"
    num_iterations_key = "num_iterations"

    # load the config file
    config_parser = ns_config.load_config_parser_from_fp(config_fp)
    score_params = ns_config.load_config_section_dict(config_parser, "score_pipeline")

    result = score_params.copy()
    if is_test:
        result[use_seed_key] = config_parser.get(test_config_section_key, use_seed_key)
        result[num_iterations_key] = config_parser.get(test_config_section_key, num_iterations_key)

    result[count_fps_or_dirs_key] = count_fps_or_dirs
    # Note: the time_prefixes_str and day_timepoints_str comma-delimited string params are not being converted to lists
    # here--that is done in the notebook, because if users run the notebooks directly, they will have to put in
    # comma-delimited strings there, so the notebooks needs to know how to deal with it.
    result[day_timepoints_str_key] = day_timepoints_str

    # the below values DO need to be converted because users have the ability to input int, float, and boolean values
    # directly into the notebooks, so the notebooks don't need to know how to convert those
    result[min_count_limit_key] = int(result[min_count_limit_key])
    result[max_fraction_acceptable_spline_density_diff_key] = float(
        result[max_fraction_acceptable_spline_density_diff_key])
    result[max_fraction_counts_excluded_key] = float(
        result[max_fraction_counts_excluded_key])
    result[use_seed_key] = bool(distutils.util.strtobool(result[use_seed_key]))
    result[num_iterations_key] = int(result[num_iterations_key])

    # if test parameters are detected, remind user! Results should not be used for real analysis
    if result[use_seed_key] or is_test:
        warnings.warn('Scoring is running in TEST MODE; do not use results for data analysis!')

    return result


def main():
    dataset_name, library_name, counts_fps_or_dirs, day_timepoints_str, is_test, config_fp = _parse_cmd_line_args()
    score_params = _set_params(counts_fps_or_dirs, day_timepoints_str, is_test, config_fp)
    full_params = ns_dcpipe.generate_notebook_params(dataset_name, library_name, score_params, config_fp)
    # Note: second argument is the *function*, not results of calling the function
    ns_pipeliner.execute_run_from_full_params(full_params, ns_dcpipe.rename_param_names_as_global_vars)


if __name__ == '__main__':
    main()