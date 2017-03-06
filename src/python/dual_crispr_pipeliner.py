# standard libraries
import argparse
import enum
import os

# ccbb libraries
from ccbbucsd.utilities.notebook_pipeliner import DATASET_NAME_KEY, ALG_NAME_KEY, execute_run

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


class DirectoryKeys(enum.Enum):
    CODE = "code"
    NOTEBOOKS = "notebooks"
    LIBRARIES = "library_defs"
    RAW_DATA = "raw"
    INTERIM_DATA = "interim"
    PROCESSED_DATA = "processed"


class PipelineSteps(enum.Enum):
    SCAFFOLD_TRIMMING = 0
    TRIMMED_READ_FILTERING = 1
    CONSTRUCT_COUNTING = 2
    COUNT_COMBINATION = 3
    COUNT_PLOTTING = 4


def parse_cmd_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_dir_name", help="name of the folder in which the fastq data to be analyzed reside")
    parser.add_argument("dataset_name", help="short human-readable name for the dataset to be analyzed")
    parser.add_argument("library_name", help="name of the construct library for the dataset to be analyzed")
    args = parser.parse_args()
    return args.fastq_dir_name, args.dataset_name, args.library_name


def generate_params(fastq_set_name, human_readable_name, num_processors, dirs_dict, spacers_file_name,
                    col_indices, min_trimmed_grna_len, max_trimmed_grna_len, full_5p_r1, full_5p_r2, full_3p_r1,
                    full_3p_r2, len_of_seq_to_match, num_allowed_mismatches):

    result = {'g_code_location': dirs_dict[DirectoryKeys.CODE],
                     'g_timestamp': "{timestamp}",
                     DATASET_NAME_KEY: human_readable_name,
                     'g_num_processors': num_processors,
                     'g_fastqs_dir': os.path.join(dirs_dict[DirectoryKeys.RAW_DATA], fastq_set_name),
                     'g_trimmed_fastqs_dir': dirs_dict["interim"],
                     'g_full_5p_r1': full_5p_r1,
                     'g_full_5p_r2': full_5p_r2,
                     'g_full_3p_r1': full_3p_r1,
                     'g_full_3p_r2': full_3p_r2,
                     'g_filtered_fastqs_dir': dirs_dict[DirectoryKeys.INTERIM_DATA],
                     'g_min_trimmed_grna_len': min_trimmed_grna_len,
                     'g_max_trimmed_grna_len': max_trimmed_grna_len,
                     'g_len_of_seq_to_match': len_of_seq_to_match,
                     ALG_NAME_KEY: "{0}mer_{1}mm_py".format(len_of_seq_to_match, num_allowed_mismatches),
                     'g_num_allowed_mismatches': num_allowed_mismatches,
                     'g_constructs_fp': os.path.join(dirs_dict[DirectoryKeys.LIBRARIES], spacers_file_name),
                     'g_col_indices_str': col_indices,
                     'g_fastq_counts_dir': dirs_dict[DirectoryKeys.PROCESSED_DATA],
                     'g_fastq_counts_run_prefix': human_readable_name + "_{timestamp}",
                     'g_collapsed_counts_dir': "{run_dir}",
                     'g_collapsed_counts_run_prefix': "{run_prefix}",
                     'g_combined_counts_dir': "{run_dir}",
                     'g_combined_counts_run_prefix': "{run_prefix}",
                     'g_plots_dir': "{run_dir}",
                     'g_plots_run_prefix': "{run_prefix}"
                     }

    return result


def run_pipeline(dirs_dict, shared_params):
    steps_to_run = [PipelineSteps.SCAFFOLD_TRIMMING,
                    PipelineSteps.TRIMMED_READ_FILTERING,
                    PipelineSteps.CONSTRUCT_COUNTING,
                    PipelineSteps.COUNT_COMBINATION,
                    PipelineSteps.COUNT_PLOTTING]

    notebooks_dir = dirs_dict[DirectoryKeys.NOTEBOOKS]
    pipeline_steps = {PipelineSteps.SCAFFOLD_TRIMMING: [notebooks_dir, "Dual CRISPR 1-Construct Scaffold Trimming.ipynb"],
                      PipelineSteps.TRIMMED_READ_FILTERING: [notebooks_dir, "Dual CRISPR 2-Constuct Filter.ipynb"],
                      PipelineSteps.CONSTRUCT_COUNTING: [notebooks_dir, "Dual CRISPR 3-Construct Counting.ipynb"],
                      PipelineSteps.COUNT_COMBINATION: [notebooks_dir, "Dual CRISPR 4-Count Combination.ipynb"],
                      PipelineSteps.COUNT_PLOTTING: [notebooks_dir, "Dual CRISPR 5-Count Plots.ipynb"]}

    execute_run(pipeline_steps, shared_params, steps_to_run, dirs_dict[DirectoryKeys.PROCESSED_DATA])

