# standard libraries
import datetime
import os

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


def get_notebooks_dir_key():
    return "notebook_dir"


def get_results_dir_key():
    return "processed_data_dir"


def get_notebook_names_list_key():
    return "notebook_basenames_list"


def get_run_prefix_key():
    return "run_prefix"


def get_run_dir_key():
    return "run_dir"


def describe_var_list(input_var_name_list):
    description_list = ["{0}: {1}\n".format(name, eval(name)) for name in input_var_name_list]
    return "".join(description_list)


def check_or_set(input_val, output_val):
    return input_val if input_val else output_val


def generate_run_prefix_and_dir_dict(parent_dir, *args, test_val=None):
    run_prefix = generate_run_prefix(*args, test_val=test_val)
    run_dir = os.path.join(parent_dir, run_prefix)
    result = {get_run_prefix_key(): run_prefix,
              get_run_dir_key(): run_dir}
    return result


def generate_run_prefix(*args, test_val=None):
    pieces = list(args)
    pieces.append(_get_timestamp(test_val))
    delimiter = _get_delimiter()
    return delimiter.join(pieces)


def generate_run_dir(results_dir, run_prefix):
    return os.path.join(results_dir, run_prefix)


def strip_run_prefix(string_to_strip, run_prefix):
    run_prefix_then_separator = run_prefix + _get_delimiter()
    result = string_to_strip.replace(run_prefix_then_separator, "")  # remove run prefix followed by separator

    separator_then_run_prefix = _get_delimiter() + run_prefix
    result = result.replace(separator_then_run_prefix, "")  # remove separator followed by run prefix

    result = result.replace(run_prefix, "")  # if run prefix exists w/o separator, remove that too

    result = result.strip()
    if not result:
        raise ValueError(
            "Stripping run prefix '{0}' from string '{1}' produces an empty string".format(run_prefix, string_to_strip))

    return result


def _get_delimiter():
    return "_"


def _get_timestamp(test_val=None):
    result = test_val  # only purpose of test_val is to make timestamps testable
    if result is None:
        result = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
    return result


