# standard libraries
import argparse

# project-specific libraries
import dual_crispr_pipeliner as ns_dcpipe

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def parse_cmd_line_args():
    # examples:
    # human_readable_name = "20160627HeLaA549"
    # library_name = "CV4"
    # day_timepoints_str = "3,14,20,28"

    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_name", help="short, alphanumeric human-readable name for the dataset to be analyzed")
    parser.add_argument("library_name", help="name of the construct library for the dataset to be analyzed")
    parser.add_argument("counts_fp_or_dir", help="the path to the file of counts to be scored, or to the directory in which the *_combined_counts.txt file to be scored resides")
    parser.add_argument("day_timepoints_str", help="comma-separated string containing the days (in order) at which data were collected")
    args = parser.parse_args()
    return args.dataset_name, args.library_name, args.counts_fp_or_dir, args.day_timepoints_str


def main():
    human_readable_name, library_name, counts_fp_or_dir, day_timepoints_str = parse_cmd_line_args()
    dirs_dict = ns_dcpipe.generate_dirs_dict()
    expt_name = ns_dcpipe.get_expt_name(human_readable_name, library_name)
    library_tuple = ns_dcpipe.get_library_params(library_name)
    params = ns_dcpipe.generate_score_params(expt_name, counts_fp_or_dir, day_timepoints_str, dirs_dict, library_tuple)

    ns_dcpipe.run_score_pipeline(dirs_dict, params)

if __name__ == '__main__':
    main()