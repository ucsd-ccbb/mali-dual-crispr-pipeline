# standard libraries
import collections
import csv
import os
import timeit

# third-party libraries
import cutadapt.scripts.cutadapt

# ccbb libraries
from ccbbucsd.utilities.basic_fastq import FastqHandler, paired_fastq_generator

__author__ = 'Amanda Birmingham'
__version__ = "0.0.1"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def run_cutadapt_global_trim(input_fastq_fp, seq_to_trim, is_5p):
    end_switch = "-g"
    end_name = "5"
    if not is_5p:
        end_switch = "-a"
        end_name = "3"

    output_fastq_fp = input_fastq_fp.replace(".fastq", "_trimmed{0}.fastq".format(end_name))
    args = [end_switch, seq_to_trim, "-o", output_fastq_fp, input_fastq_fp]
    cutadapt.scripts.cutadapt.main(args)
    return output_fastq_fp


def trim_linked_scaffold(fastq_fp, scaffold_seq_5p, scaffold_seq_3p):
    output_fastq_fp = fastq_fp.replace(".fastq", "_trimmed53.fastq")
    args = ["-a", "{0}...{1}".format(scaffold_seq_5p,scaffold_seq_3p), "-o", output_fastq_fp, fastq_fp]
    cutadapt.scripts.cutadapt.main(args)
    return output_fastq_fp


def trim_scaffold(fastq_fp, scaffold_seq_5p=None, scaffold_seq_3p=None):
    curr_fastq_fp = fastq_fp

    if scaffold_seq_5p is not None:
        curr_fastq_fp = run_cutadapt_global_trim(curr_fastq_fp, scaffold_seq_5p, True)

    if scaffold_seq_3p is not None:
        curr_fastq_fp = run_cutadapt_global_trim(curr_fastq_fp, scaffold_seq_3p, False)

    return curr_fastq_fp


def test_cutadapt(fw_fastq_fp, rv_fastq_fp, output_dir):
    dongxin_5p_r1 = "GTGGAAAGGACGAAACACCG"
    dongxin_5p_r2 = "GCTATTTCTAGCTCTAAAAC"
    full_5p_r1 = "TATATATCTTGTGGAAAGGACGAAACACCG"
    full_5p_r2 = "CCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC"
    full_3p_r1 = "GTTTCAGAGCTATGCTGGAAACTGC"
    full_3p_r2 = "CAAACAAGGCTTTTCTCCAAGG"

    # fw_trimmed_fastq_fp = trim_scaffold(fw_fastq_fp, full_5p_r1, full_3p_r1)
    # rv_trimmed_fastq_fp = trim_scaffold(rv_fastq_fp, full_5p_r2, full_3p_r2)
    # test_name = "full_5p_full_3p_2_global_trims_Hela-CV4-d3-1_S1_L001"

    # fw_trimmed_fastq_fp = trim_scaffold(fw_fastq_fp, dongxin_5p_r1, full_3p_r1)
    # rv_trimmed_fastq_fp = trim_scaffold(rv_fastq_fp, dongxin_5p_r2, full_3p_r2)
    # test_name = "dongxin_5p_full_3p_2_global_trims_Hela-CV4-d3-1_S1_L001"

    # fw_trimmed_fastq_fp = trim_linked_scaffold(fw_fastq_fp, full_5p_r1, full_3p_r1)
    # rv_trimmed_fastq_fp = trim_linked_scaffold(rv_fastq_fp, full_5p_r2, full_3p_r2)
    # test_name = "full_5p_full_3p_2_linked_trims_Hela-CV4-d3-1_S1_L001"

    # fw_trimmed_fastq_fp = trim_scaffold(fw_fastq_fp, full_5p_r1, full_3p_r1)
    # rv_trimmed_fastq_fp = trim_scaffold(rv_fastq_fp, full_5p_r2, full_3p_r2)
    # test_name = "full_5p_full_3p_2_global_trims_Hela-CV4-d14-1_S3_L001"

    # fw_trimmed_fastq_fp = trim_scaffold(fw_fastq_fp, dongxin_5p_r1, full_3p_r1)
    # rv_trimmed_fastq_fp = trim_scaffold(rv_fastq_fp, dongxin_5p_r2, full_3p_r2)
    # test_name = "dongxin_5p_full_3p_2_global_trims_Hela-CV4-d14-1_S3_L001"

    fw_trimmed_fastq_fp = trim_linked_scaffold(fw_fastq_fp, full_5p_r1, full_3p_r1)
    rv_trimmed_fastq_fp = trim_linked_scaffold(rv_fastq_fp, full_5p_r2, full_3p_r2)
    test_name = "full_5p_full_3p_linked_trims_Hela-CV4-d14-1_S3_L001"


    output_fp = os.path.join(output_dir, "{0}.csv".format(test_name))
    counters_dict = examine_cutadapt_results(fw_trimmed_fastq_fp, rv_trimmed_fastq_fp)
    for curr_key, curr_item in counters_dict.items():
        write_counter(curr_key, curr_item, output_fp, test_name)

    return test_name, output_fp


def examine_cutadapt_results(fw_trimmed_fp, rv_trimmed_fp):
    result = {}
    fw_fastq_handler = FastqHandler(fw_trimmed_fp)
    rv_fastq_handler = FastqHandler(rv_trimmed_fp)

    num_fastq_pairs = 0
    r1_counter = collections.Counter()
    r2_counter = collections.Counter()
    shared_counter = collections.Counter()

    paired_fastq_seqs = paired_fastq_generator(fw_fastq_handler, rv_fastq_handler)
    for curr_pair_seqs in paired_fastq_seqs:
        num_fastq_pairs += 1
        if num_fastq_pairs % 10000 == 0:
            print("On fastq pair number {0}".format(num_fastq_pairs))
            if num_fastq_pairs == 300000:
                print("Ending")
                break

        r1_len = len(curr_pair_seqs[0])
        r2_len = len(curr_pair_seqs[1])
        r1_counter[r1_len] += 1
        r2_counter[r2_len] += 1

        if r1_len == r2_len == 20:
            shared_counter["perfect"] += 1
        elif (r1_len == 19 and r2_len == 20) or (r1_len == 20 and r1_len == 19):
            shared_counter["close"] += 1

    result["r1"] = r1_counter
    result["r2"] = r2_counter
    result["shared"] = shared_counter
    return result


def write_counter(counter_header, counter, output_fp, test_name):
    with open(output_fp, 'a') as file_handle:
        writer = csv.writer(file_handle)
        writer.writerow(["#", "{0}_{1}".format(counter_header, test_name)])
        for key, value in counter.items():
            writer.writerow([key, value])


def cutadapt_tester():
    project_dir = "/Users/Birmingham/Repositories/ccbb_tickets/20160210_mali_crispr/data"

    # logging.basicConfig(filename=os.path.join(project_dir,'construct_counting2_tester11.log'), level=logging.DEBUG)

    # grnas_fp = os.path.join(project_dir, "raw/grna_name_by_seq.txt")
    # constructs_fp = os.path.join(project_dir, "raw/CV4_2spacers.txt")
    fw_fastq_fp = os.path.join(project_dir, "raw/20160403_data/Hela-CV4-d14-1-34707397/Hela-CV4-d14-1_S3_L001_R1_001.fastq")
    rv_fastq_fp = os.path.join(project_dir, "raw/20160403_data/Hela-CV4-d14-1-34707397/Hela-CV4-d14-1_S3_L001_R2_001.fastq")
    output_dir = os.path.join(project_dir, "processed")

    start_time = timeit.default_timer()
    print("start time: {0}".format(start_time))

    # fw_fastq_fp = os.path.join(project_dir, "raw/20160506_data/U2OS-CAS9TREX2A-t14_S1_L001_R1_001.fastq")
    # rv_fastq_fp = os.path.join(project_dir, "raw/20160506_data/U2OS-CAS9TREX2A-t14_S1_L001_R2_001.fastq")
    test_name, output_fp = test_cutadapt(fw_fastq_fp, rv_fastq_fp, output_dir)

    end_time = timeit.default_timer()
    elapsed_time = end_time - start_time
    time_msg = "{0}: {1}".format(test_name, elapsed_time)
    print(time_msg)
    with open(output_fp, 'a') as file_handle:
        file_handle.write(time_msg)

cutadapt_tester()
