import configparser


def load_config_parser_from_fp(config_fp):
    # Note: configparser has a read_file method that could be used instead, but  structuring the code
    # this way makes it easier to unit-test bc can pass test string to load_config_settings
    with open(config_fp, 'r') as config_file:
        config_string = config_file.read()

    return load_config_parser(config_string)


def load_config_parser(config_string):
    result = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    result.read_string(config_string)
    return result


def load_config_section_dict(config_parser, section_name=None):
    result = {}

    if section_name is None:
        section_name = config_parser.default_section

    for curr_key, _ in config_parser.items(section_name):
        result[curr_key] = config_parser.get(section_name, curr_key)

    return result


# # May need this someday for windows support?
# def convert_param_strings(config_params):
#     import os
#     for curr_key, curr_value in config_params.items():
#         if curr_key.endswith("_dir"):
#             config_params[curr_key] = os.path.normpath(curr_value)
#     return config_params
