# standard libraries
import warnings

# ccbb libraries
import ccbbucsd.utilities.files_and_paths as ns_file

_COMMENT_CHAR = "#"
_LIBRARY_NAME_KEY = "library_name"
_LIBRARY_FP_KEY = "library_fp"
_LIBRARY_SETTINGS_KEYS = ["max_trimmed_grna_len", "min_trimmed_grna_len"]
_SETTINGS_SEPARATOR = "="


def id_library_info(library_name, library_dir):
    libraries_dict_list = _extract_library_info_from_files(library_name, library_dir)
    result = _id_and_validate_library(library_name, libraries_dict_list)
    return result


def _extract_library_info_from_files(library_name, library_dir):
    libraries_dict_list = []

    text_fps = ns_file.get_filepaths_from_wildcard(library_dir, ".txt")
    for curr_fp in text_fps:
        with open(curr_fp, 'r') as curr_file:
            library_dict = _mine_library_file(library_name, curr_fp, curr_file)
            if library_dict is not None:
                libraries_dict_list.append(library_dict)

    return libraries_dict_list


def _mine_library_file(library_name, curr_fp, curr_file):
    result = None

    first_line = curr_file.readline()
    first_line_settings = _id_and_trim_settings_line(first_line)
    if first_line_settings is not None:
        if first_line_settings[0] == _LIBRARY_NAME_KEY:
            if first_line_settings[1] == library_name:
                result = {first_line_settings[0]: first_line_settings[1],
                          _LIBRARY_FP_KEY: curr_fp}

                # read as many lines as there should be settings
                for i in range(0, len(_LIBRARY_SETTINGS_KEYS)):
                    curr_line = curr_file.readline()
                    curr_settings = _id_and_trim_settings_line(curr_line)
                    if curr_settings is not None:
                        result[curr_settings[0]] = curr_settings[1]

    return result


def _id_and_trim_settings_line(a_line):
    result = None  # assume line does not contain valid settings

    if a_line.startswith(_COMMENT_CHAR):
        # Take any comment characters off the left end of the string, and then whitespace off both ends
        trimmed_line = a_line.lstrip(_COMMENT_CHAR).strip()
        # split on = and trim result
        split_line = [x.strip() for x in trimmed_line.split(_SETTINGS_SEPARATOR)]

        # a real settings line should split into two pieces (key and setting) on the settings separator,
        # and neither of the two pieces should be empty strings.  Note that in python an empty string
        # evaluates to false, so all(split_line) "ands" together every string in split_line, and we are
        # checking that the "and" of all these strings is true--i.e., every string is non-empty
        if len(split_line) == 2 and all(split_line):
            result = tuple(split_line)

    return result


def _id_and_validate_library(library_name, libraries_dict_list):
    # validate that we found one and only one library by this name
    _validate_library_id(library_name, libraries_dict_list)

    library_dict = libraries_dict_list[0]

    # validate that all expected keys are present
    _validate_settings_keys(library_name, library_dict)

    # validate that all keys have values
    _validate_settings_values(library_name, library_dict)

    return library_dict


def _validate_library_id(library_name, libraries_dict_list):
    # validate that we found one and only one library by this name
    if len(libraries_dict_list) == 0:
        raise ValueError("No library file found with library name '{0}'".format(library_name))
    elif len(libraries_dict_list) > 1:
        raise ValueError("Multiple library files found with library name '{0}': {1}".format(
            library_name, ", ".join([x[_LIBRARY_FP_KEY] for x in libraries_dict_list])))


def _validate_settings_keys(library_name, library_dict):
    expected_keys_list = _LIBRARY_SETTINGS_KEYS
    expected_keys_list.extend([_LIBRARY_NAME_KEY, _LIBRARY_FP_KEY])
    expected_keys_set = set(expected_keys_list)
    found_keys_set = set(library_dict.keys())

    # warnings done first so that we see them even if there are *also* errors later
    ignored_keys = found_keys_set - expected_keys_set
    if len(ignored_keys) != 0:
        warnings.warn("Library '{0}' includes the following ignored settings: {1}".format(
            library_name, ", ".join(ignored_keys)))

    missing_keys = expected_keys_set - found_keys_set
    if len(missing_keys) != 0:
        raise ValueError("Library '{0}' is missing the following expected settings: {1}".format(
            library_name, ", ".join(missing_keys)))


def _validate_settings_values(library_name, library_dict):
    problem_keys = []
    for curr_key, curr_value in library_dict.items():
        if not curr_value:
            problem_keys.append(curr_key)

    if len(problem_keys) > 0:
        raise ValueError("Library '{0}' is missing settings values for the following settings: {1}".format(
            library_name, ", ".join(problem_keys)))


