# standard libraries
from os import path

# third-party libraries
import pandas

# project-specific libraries
from ccbbucsd.malicrispr.experiment_constants import ExperimentConstants
from ccbbucsd.malicrispr.dual_crispr_headers_manager import DualCrisprHeadersManager
from ccbbucsd.malicrispr.dual_crispr_file_manager import DualCrisprFileManager
from ccbbucsd.malicrispr.median_normalization import add_normalization_across_expts
from ccbbucsd.malicrispr.fold_change import add_timept_and_plasmid_fold_changes

__author__ = 'Amanda Birmingham'
__version__ = "0.2.0"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def init_constants():
    project_path = "/Users/Birmingham/Repositories/custom_analysis_pipelines/20160210_mali_crispr/data/"

    # Mali test dataset 1
    dataset_name = "mali_test_dataset1"
    working_dir = path.join(project_path, "test_data/dataset1")
    raw_counts_fp = path.join(working_dir, "dataset1_raw_counts.txt")
    plasmid_fp = None
    count_headers = ["expt1_day3", "expt1_day16"]
    experiment_set_prefixes = ["expt1"]
    expected_num_constructs = 22910
    earliest_timept = "day3"
    negative_control_genes = ['NonTargetingControlGuideForHuman']

    return ExperimentConstants(dataset_name, working_dir, raw_counts_fp, plasmid_fp, count_headers,
                               experiment_set_prefixes, earliest_timept, negative_control_genes,
                               expected_num_constructs)


def normalize_and_save(file_manager, raw_counts_fp, plasmid_fp, num_expected_constructs):
    working_df = file_manager.get_counts_dataframe(raw_counts_fp, plasmid_fp)
    add_normalization_across_expts(working_df, file_manager.headers_mgr, num_expected_constructs)
    norm_fp = file_manager.write_normalized_to_file(working_df)
    return norm_fp


def calc_fold_changes_and_save(file_manager, df_csv_fp, expt_set_prefixes, earliest_timept):
    working_df = pandas.read_csv(df_csv_fp)
    add_timept_and_plasmid_fold_changes(working_df, file_manager, expt_set_prefixes, earliest_timept)
    fc_fp = file_manager.write_foldchanges_to_file(working_df)
    return fc_fp


def calc_fitness_scores_and_save(file_manager, df_csv_fp):
    working_df = pandas.read_csv(df_csv_fp)
    # TODO: Pick up here!!
    # add_fitness_scores_and_p_values()
    # fitness_fp = file_manager.write_fitnesses_to_file(working_df)
    # return fitness_fp


def run_pipeline(expt_constants):
    headers_mgr = DualCrisprHeadersManager()
    file_manager = DualCrisprFileManager(headers_mgr, expt_constants)
    norm_csv_fp = normalize_and_save(file_manager, expt_constants.raw_counts_fp, expt_constants.plasmid_fp,
                                     expt_constants.expected_num_constructs)
    fc_fp = calc_fold_changes_and_save(file_manager, norm_csv_fp, expt_constants.experiment_set_prefixes,
                                       expt_constants.earliest_timept)
    # fitnesses_fp = calc_fitness_scores_and_save(file_manager, fc_fp)


g_constants = init_constants()
run_pipeline(g_constants)
