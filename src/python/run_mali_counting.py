# standard libraries
import argparse

# project-specific libraries
import dual_crispr_pipeliner as ns_dcpipe

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def parse_cmd_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_dir_name", help="name of the folder in which the fastq data to be analyzed reside")
    parser.add_argument("dataset_name", help="short, alphanumeric human-readable name for the dataset to be analyzed")
    parser.add_argument("library_name", help="name of the construct library for the dataset to be analyzed")
    args = parser.parse_args()
    return args.fastq_dir_name, args.dataset_name, args.library_name


def main():
    fastq_set_name, human_readable_name, library_name = parse_cmd_line_args()
    dirs_dict = ns_dcpipe.generate_dirs_dict(fastq_set_name=fastq_set_name)
    expt_name = ns_dcpipe.get_expt_name(human_readable_name, library_name)
    library_tuple = ns_dcpipe.get_library_params(library_name)
    count_params = ns_dcpipe.generate_counts_params(fastq_set_name, expt_name, dirs_dict, library_tuple)

    ns_dcpipe.run_counts_pipeline(dirs_dict, count_params)

if __name__ == '__main__':
    main()