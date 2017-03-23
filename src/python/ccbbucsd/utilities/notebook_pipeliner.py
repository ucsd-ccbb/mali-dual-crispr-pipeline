# standard libraries
import os
import re

# ccbb libraries
import ccbbucsd.utilities.analysis_run_prefixes as ns_runs
import ccbbucsd.utilities.notebook_runner as ns_notebook

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


def execute_run_from_full_params(run_params, params_preprocess_func=None, test_val=None):
    notebook_dir = run_params[ns_runs.get_notebooks_dir_key()]
    notebook_basenames_list = run_params[ns_runs.get_notebook_names_list_key()]
    results_dir = run_params[ns_runs.get_results_dir_key()]

    _execute_run(notebook_dir, notebook_basenames_list, results_dir, run_params, params_preprocess_func, test_val)


def _execute_run(notebook_dir, notebook_basenames_list, results_dir, run_params, params_preprocess_func=None,
                 test_val=None):

    run_prefix, run_dir_name, extended_params = _add_run_prefix_and_dir_to_params(results_dir, run_params,
                                                                                  params_preprocess_func, test_val)

    methods_dir = _create_run_and_methods_dirs(run_dir_name)

    for curr_notebook_basename in notebook_basenames_list:
        _run_and_output_notebook(notebook_dir, curr_notebook_basename, extended_params, run_prefix, methods_dir)


def _add_run_prefix_and_dir_to_params(results_dir, run_params, params_preprocess_func=None, test_val=None):
    if not ns_runs.get_run_prefix_key() in run_params:
        run_params[ns_runs.get_run_prefix_key()] = ns_runs.generate_run_prefix(test_val=test_val)
    run_prefix = run_params[ns_runs.get_run_prefix_key()]

    if not ns_runs.get_run_dir_key() in run_params:
        run_params[ns_runs.get_run_dir_key()] = ns_runs.generate_run_dir(results_dir, run_prefix)
    run_dir = run_params[ns_runs.get_run_dir_key()]

    if params_preprocess_func is not None:
        run_params = params_preprocess_func(run_params)

    return run_prefix, run_dir, run_params


def _run_and_output_notebook(notebook_dir, base_notebook_filename, params_dict, run_prefix, methods_dir):
    notebook_out_fp = _get_output_fp(base_notebook_filename, run_prefix, methods_dir, ".ipynb")
    ns_notebook.execute_notebook(base_notebook_filename, params_dict, notebook_out_fp, run_path=notebook_dir)
    html_out_fp = ns_notebook.export_notebook_to_html(notebook_out_fp, methods_dir)
    return notebook_out_fp, html_out_fp


def _create_run_and_methods_dirs(run_dir_name):
    methods_dir_name = os.path.join(run_dir_name, _get_methods_folder_name())
    # makedirs differs from mkdir in that it will make intermediate directories that don't already exist--
    # in this case, it will make the run dir that is the parent of the methods dir
    os.makedirs(methods_dir_name, exist_ok=True)  # True = is OK if path already exists
    return methods_dir_name


def _get_methods_folder_name():
    return "methods"


def _get_output_fp(notebook_fp, run_prefix, methods_dir, output_ext):
    _, notebook_name = os.path.split(notebook_fp)
    new_base = _mangle_notebook_name(notebook_name, run_prefix)
    new_fp = os.path.join(methods_dir, new_base + output_ext)
    return new_fp


def _mangle_notebook_name(notebook_filename, run_prefix):
    delimiter = "_"
    name_base, _ = os.path.splitext(notebook_filename)
    lower_base = name_base.lower().strip()
    delimited_notebook_base = re.sub("\s+", delimiter, lower_base)
    new_base = "{0}{1}{2}".format(run_prefix, delimiter, delimited_notebook_base)
    return new_base
