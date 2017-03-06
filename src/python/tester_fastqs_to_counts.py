# standard libraries
import logging
import os
import timeit

# project-specific libraries
from ccbbucsd.malicrispr.construct_counter import generate_construct_counts

__author__ = 'Amanda Birmingham'
__version__ = "0.0.1"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def tester():
    project_dir = "/Users/Birmingham/Repositories/ccbb_tickets/20160210_mali_crispr/data"

    #logging.basicConfig(filename=os.path.join(project_dir,'construct_counting2_tester11.log'), level=logging.DEBUG)

    grnas_fp = os.path.join(project_dir, "raw/grna_name_by_seq.txt")
    constructs_fp = os.path.join(project_dir, "raw/CV4_2spacers.txt")
    fw_fastq_fp = os.path.join(project_dir, "raw/20160403_data/Hela-CV4-d14-1-34707397/full_5p_full_3p_linked_trims/Hela-CV4-d14-1_S3_L001_R1_001_trimmed53.fastq")
    rv_fastq_fp = os.path.join(project_dir, "raw/20160403_data/Hela-CV4-d14-1-34707397/full_5p_full_3p_linked_trims/Hela-CV4-d14-1_S3_L001_R2_001_trimmed53.fastq")
    # fw_fastq_fp = os.path.join(project_dir, "raw/20160506_data/U2OS-CAS9TREX2A-t14_S1_L001_R1_001.fastq")
    # rv_fastq_fp = os.path.join(project_dir, "raw/20160506_data/U2OS-CAS9TREX2A-t14_S1_L001_R2_001.fastq")
    logfile_fp = os.path.join(project_dir, "processed/Hela-CV4-d14-1_S3_L001_trimmed53_counting.log")
    output_fp = os.path.join(project_dir, "processed/Hela-CV4-d14-1_S3_L001_trimmed53.txt")

    generate_construct_counts(grnas_fp, constructs_fp, logfile_fp, output_fp, fw_fastq_fp, rv_fastq_fp)

    # start_time = timeit.default_timer()
    # counts_info_tuple = match_and_count_constructs(grnas_fp, constructs_fp, fw_fastq_fp, rv_fastq_fp)
    # end_time = timeit.default_timer()
    # elapsed_time = end_time - start_time
    # print("0mm anywhere matcher CV4-Plasmid_S1_L001_R2_001 elapsed time: {0}".format(elapsed_time))
    # counts_by_construct = counts_info_tuple[0]
    # counts_summary_list = counts_info_tuple[1]
    # write_counts(counts_by_construct, counts_summary_list, output_fp)


tester()

import enum


# class LoggerTypes(enum.Enum):
#     INELIGIBLE_READ = 1
#     UNRECOGNIZED_READ = 2
#     UNRECOGNIZED_CONSTRUCT = 3
#
#
# def _set_up_loggers(output_dir, run_prefix):
#     for logger_type in LoggerTypes:
#         logger_name = logger_type.name
#         curr_logger = logging.getLogger(logger_type)
#
#         log_fp = os.path.join(output_dir, "{0}_{1}.log".format(run_prefix, logger_name))
#         fh = logging.FileHandler(log_fp)
#         formatter = logging.Formatter('%(message)s')
#         fh.setFormatter(formatter)
#         curr_logger.addHandler(fh)