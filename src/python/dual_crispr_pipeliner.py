# standard libraries
import enum
import os
import warnings

# ccbb libraries
import ccbbucsd.utilities.files_and_paths as ns_file
from ccbbucsd.utilities.notebook_pipeliner import DATASET_NAME_KEY, ALG_NAME_KEY, execute_run

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


_NUM_PROCESSORS = 7 # assumes on an 8-vCPU instance


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
    SCORING_PREP = 5
    ABUNDANCE_THRESHOLD_PICKING = 6
    SCORING = 7


# TODO: either refactor out directory-making or rename method to alert users to side-effects
def generate_dirs_dict(main_dir="/home/ec2-user/mali-dual-crispr-pipeline", data_dir="/data", fastq_set_name=""):
    dirs_dict = {}
    dirs_dict[DirectoryKeys.CODE] = os.path.join(main_dir, "src", "python")
    dirs_dict[DirectoryKeys.NOTEBOOKS] = os.path.join(main_dir, DirectoryKeys.NOTEBOOKS.value)
    dirs_dict[DirectoryKeys.LIBRARIES] =  os.path.join(main_dir, "library_definitions")
    dirs_dict[DirectoryKeys.RAW_DATA] = os.path.join(data_dir, DirectoryKeys.RAW_DATA.value)
    dirs_dict[DirectoryKeys.INTERIM_DATA] = os.path.join(data_dir, DirectoryKeys.INTERIM_DATA.value, fastq_set_name)
    dirs_dict[DirectoryKeys.PROCESSED_DATA] = os.path.join(data_dir, DirectoryKeys.PROCESSED_DATA.value)

    ns_file.verify_or_make_dir(dirs_dict[DirectoryKeys.RAW_DATA])
    ns_file.verify_or_make_dir(dirs_dict[DirectoryKeys.INTERIM_DATA])
    ns_file.verify_or_make_dir(dirs_dict[DirectoryKeys.PROCESSED_DATA])

    return dirs_dict


def get_library_params(library_name):
    # you'll change these lines only when changing the construct library from which you're analyzing data.
    # Note that the file named here is expected to be in tab-delimited text format with no quotes around text values,
    # and that the first line will be skipped (as a header row).  For this file, provide the zero-based indexes of the
    # 7 columns for the following values, in the order given here:
    # Construct ID (probe pair id), e.g. FGFR3_chr4_1795735__SMAD4_chr18_48586241
    # Gene A name (target A name), e.g. FGFR3
    # Probe A name, e.g. FGFR3_chr4_1795735
    # Probe A sequence, e.g. GACGCGCTGCTCCGTCCCCA
    # Gene B name (target B name), e.g. SMAD4
    # Probe B name, e.g. SMAD4_chr18_48586241
    # Probe B sequence, e.g. AATGCAAGCTCATTGTGAAC

    if library_name == "CV4":
        # for CV4 library:
        spacers_file_name = "CV4_2spacers_w_probe_names_wo_duplicate.txt"
        col_indices = "1,3,6,7,8,11,12"
        min_trimmed_grna_len = 19
        max_trimmed_grna_len = 21
    else:
        raise ValueError("Unrecognized library name: {0}".format(library_name))

    # for MV4 library:
    # spacers_file_name = "Metabolism_dual_spacers.txt"
    # col_indices = "1,6,10"
    # min_trimmed_grna_len = 19
    # max_trimmed_grna_len = 21

    # for CRV4 library:
    # spacers_file_name = "Cancer_rep2_CRV4.txt"
    # col_indices = "1,3,5"
    # min_trimmed_grna_len = 19
    # max_trimmed_grna_len = 24

    return spacers_file_name, col_indices, min_trimmed_grna_len, max_trimmed_grna_len, library_name


def get_expt_name(human_readable_name, library_name):
    if not human_readable_name.isalnum():
        raise ValueError("Human-readable name '{0}' is invalid; must contain only letters and numbers.".format(
            human_readable_name))

    if not library_name.isalnum():
        raise ValueError("Library name '{0}' is invalid; must contain only letters and numbers.".format(
            library_name))

    return "{0}{1}".format(human_readable_name, library_name)


def generate_counts_params(fastq_set_name, expt_name, dirs_dict, library_tuple):
    # *********************************************************
    # you'll change these lines only when changing the scaffold sequence with which your library was constructed.
    # Note that these should be the FULL scaffold sequence, even if the sequencing will not capture the whole 3'
    # scaffold sequence; the cutadapt settings used will still identify a significantly truncated 3' match.
    full_5p_r1 = "TATATATCTTGTGGAAAGGACGAAACACCG"
    full_5p_r2 = "CCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC"
    full_3p_r1 = "GTTTCAGAGCTATGCTGGAAACTGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTACTGAG"
    full_3p_r2 = "CAAACAAGGCTTTTCTCCAAGGGATATTTATAGTCTCAAAACACACAATTACTTTACAGTTAGGGTGAGTTTCCTTTTGTGCTGTTTTTTAAAATA"

    # *********************************************************
    # you'll probably decide on these values before beginning screen analysis, set them once, and then leave them alone
    # unless you change your analysis approach and/or hardware
    len_of_seq_to_match = 19
    num_allowed_mismatches = 1

    spacers_file_name = library_tuple[0]
    col_indices = library_tuple[1]
    min_trimmed_grna_len = library_tuple[2]
    max_trimmed_grna_len = library_tuple[3]


    result = {'g_code_location': dirs_dict[DirectoryKeys.CODE],
                     'g_timestamp': "{timestamp}",
                     DATASET_NAME_KEY: expt_name,
                     'g_num_processors': _NUM_PROCESSORS,
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
                     'g_fastq_counts_run_prefix': expt_name + "_{timestamp}",
                     'g_collapsed_counts_dir': "{run_dir}",
                     'g_collapsed_counts_run_prefix': "{run_prefix}",
                     'g_combined_counts_dir': "{run_dir}",
                     'g_combined_counts_run_prefix': "{run_prefix}",
                     'g_plots_dir': "{run_dir}",
                     'g_plots_run_prefix': "{run_prefix}"
                     }

    return result


def run_counts_pipeline(dirs_dict, shared_params):
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


def generate_score_params(expt_name, counts_fp_or_dir, day_timepoints_str, dirs_dict, library_tuple, is_test):
    spacers_file_name = library_tuple[0]
    col_indices = library_tuple[1]

    use_seed = False
    num_iterations = 1000
    if is_test:
        use_seed = True
        num_iterations = 2
        warnings.warn('Scoring is running in TEST MODE; do not use results for data analysis!')


    result = {'g_code_location': dirs_dict[DirectoryKeys.CODE],
                 'g_timestamp': "{timestamp}",
                 DATASET_NAME_KEY: expt_name,
                 'g_num_processors': _NUM_PROCESSORS,
                 ALG_NAME_KEY: "method2",
                 'g_constructs_fp': os.path.join(dirs_dict[DirectoryKeys.LIBRARIES], spacers_file_name),
                 'g_col_indices_str': col_indices,
                 'g_count_fps_or_dirs': counts_fp_or_dir,
                 'g_prepped_counts_run_prefix': expt_name + "_{timestamp}",
                 'g_prepped_counts_dir': "{run_dir}",
                 'g_min_count_limit': 10,  # Note: in absolute counts, not log2.
                 'g_max_fraction_acceptable_spline_density_diff': 0.02,  # % of diff between max spline and min density
                 'g_max_fraction_counts_excluded': 0.95,  # any threshold throwing out >x% of counts is not acceptable
                 'g_thresholds_dir': "{run_dir}",
                 'g_thresholds_run_prefix': "{run_prefix}",
                 'g_collapsed_counts_run_prefix': "{run_prefix}",
                 'g_use_seed': use_seed,  # This should be *False* when running on real data!
                 'g_num_iterations': num_iterations,  # This should be *1000* when running on real data!
                 'g_counts_fp': "{run_dir}",
                 'g_day_timepoints_str': day_timepoints_str,
                 'g_scoring_dir': "{run_dir}"
          }

    return result


def run_score_pipeline(dirs_dict, shared_params):
    steps_to_run = [PipelineSteps.SCORING_PREP,
                    PipelineSteps.ABUNDANCE_THRESHOLD_PICKING,
                    PipelineSteps.SCORING]

    notebooks_dir = dirs_dict[DirectoryKeys.NOTEBOOKS]
    pipeline_steps = {PipelineSteps.SCORING_PREP: [notebooks_dir, "Dual CRISPR 6-Scoring Preparation.ipynb"],
                      PipelineSteps.ABUNDANCE_THRESHOLD_PICKING: [notebooks_dir,
                                                                  "Dual CRISPR 7-Abundance Thresholds.ipynb"],
                      PipelineSteps.SCORING: [notebooks_dir, "Dual CRISPR 8-Construct Scoring.ipynb"]}

    execute_run(pipeline_steps, shared_params, steps_to_run, dirs_dict[DirectoryKeys.PROCESSED_DATA])

