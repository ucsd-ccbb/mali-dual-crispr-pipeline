# standard libraries
import distutils.util
import enum
import os

# ccbb libraries
import ccbb_pyutils.analysis_run_prefixes as ns_runs
import ccbb_pyutils.config_loader as ns_config
import ccbb_pyutils.files_and_paths as ns_file
import ccbb_pyutils.notebook_runner as ns_notebook
import ccbb_pyutils.string_utils as ns_strings

# project-specific libraries
import dual_crispr.constructs_settings_extracter as ns_settings

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


class DirectoryKeys(enum.Enum):
    NOTEBOOKS = "notebook_dir"
    LIBRARIES = "libraries_dir"
    RAW_DATA = "raw_data_dir"
    INTERIM_DATA = "interim_data_dir"
    PROCESSED_DATA = "processed_data_dir"


def get_config_fp_or_default(config_fp):
    if config_fp is None:
        config_fp = os.path.join(os.environ['HOME'], "dual_crispr", "config.txt")
    return config_fp


def set_up_data_subdirs_and_get_machine_configs(config_fp=None):
    config_fp = get_config_fp_or_default(config_fp)
    result = get_machine_config_params(config_fp)
    _verify_or_make_data_subdirs(result)
    return result


def get_machine_config_params(config_fp=None):
    keep_gz_key = "keep_gzs"
    num_processors_key = "num_processors"

    config_fp = get_config_fp_or_default(config_fp)
    config_parser = ns_config.load_config_parser_from_fp(config_fp)

    machine_config_key = config_parser.get(config_parser.default_section, "machine_configuration")
    machine_config_dict = ns_config.load_config_section_dict(config_parser, machine_config_key)
    machine_config_dict[keep_gz_key] = bool(distutils.util.strtobool(machine_config_dict[keep_gz_key]))
    machine_config_dict[num_processors_key] = int(machine_config_dict[num_processors_key])
    return machine_config_dict


def generate_notebook_params(expt_name, library_name, arg_based_params_dict, config_fp=None):
    _validate_expt_name(expt_name)
    result = {"dataset_name": expt_name}

    config_fp = get_config_fp_or_default(config_fp)
    machine_config_dict = set_up_data_subdirs_and_get_machine_configs(config_fp)
    result.update(machine_config_dict)

    library_dict = ns_settings.id_library_info(library_name, result[DirectoryKeys.LIBRARIES.value])
    library_dict["min_trimmed_grna_len"] = int(library_dict["min_trimmed_grna_len"])
    library_dict["max_trimmed_grna_len"] = int(library_dict["max_trimmed_grna_len"])
    result.update(library_dict)

    # create and add run prefix, run directory
    run_prefix_and_dir_dict = ns_runs.generate_run_prefix_and_dir_dict(
        machine_config_dict[DirectoryKeys.PROCESSED_DATA.value], expt_name)
    result.update(run_prefix_and_dir_dict)

    # Note: the reason arg_based_params_dict is added *last* even though it could have been used as the
    # base of the params dict is that later additions can overwrite earlier ones, and the arguments params
    # should have the final say :)
    result.update(arg_based_params_dict)

    notebook_basenames_str = result[ns_runs.get_notebook_names_list_key()]
    notebook_basenames_list = ns_strings.split_delimited_string_to_list(notebook_basenames_str)
    result[ns_runs.get_notebook_names_list_key()] = notebook_basenames_list

    for curr_notebook_name in notebook_basenames_list:
        # load the notebook and extract its params
        curr_fp = os.path.join(result[DirectoryKeys.NOTEBOOKS.value], curr_notebook_name)
        notebook_params = ns_notebook.get_params_from_nb_file(curr_fp)

        for curr_param in notebook_params:
            curr_param_name = curr_param.name
            # for any param that isn't in the dict I've built up
            if curr_param_name not in result:
                if curr_param_name.endswith("_dir"):
                    result[curr_param_name] = run_prefix_and_dir_dict[ns_runs.get_run_dir_key()]
                elif curr_param_name.endswith("_run_prefix"):
                    result[curr_param_name] = run_prefix_and_dir_dict[ns_runs.get_run_prefix_key()]
                # if not one of these two, call it an armadillo and leave it alone :) (ref Just So Stories)

    return result


def rename_param_names_as_global_vars(params_dict):
    prefix = "g_"
    result = {}
    for curr_key, curr_val in params_dict.items():
        revised_key = curr_key if curr_key.startswith(prefix) else prefix + curr_key
        result[revised_key] = curr_val
    return result


def _verify_or_make_data_subdirs(machine_config_params):
    ns_file.verify_or_make_dir(machine_config_params[DirectoryKeys.RAW_DATA.value])
    ns_file.verify_or_make_dir(machine_config_params[DirectoryKeys.INTERIM_DATA.value])
    ns_file.verify_or_make_dir(machine_config_params[DirectoryKeys.PROCESSED_DATA.value])


def _validate_expt_name(human_readable_name):
    if not human_readable_name.isalnum():
        raise ValueError("Human-readable name '{0}' is invalid; must contain only letters and numbers.".format(
            human_readable_name))

