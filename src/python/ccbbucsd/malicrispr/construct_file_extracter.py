# third-party libraries
import pandas

# ccbb libraries
from ccbbucsd.utilities.bio_seq_utilities import trim_seq

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"

_CONSTRUCT_ID = "construct_id"
_PROBE_A_SEQ = "probe_a_seq"
_PROBE_B_SEQ = "probe_b_seq"
_PROBE_A_NAME = "probe_a_id"
_PROBE_B_NAME = "probe_b_id"
_TARGET_A_NAME = "target_a_id"
_TARGET_B_NAME = "target_b_id"
_TARGET_PAIR_ID = "target_pair_id"
_PROBE_PAIR_ID = "probe_pair_id"
_HEADER_DIVIDER = "_"


def get_potential_annotation_headers():
    return [_CONSTRUCT_ID, _PROBE_A_SEQ, _PROBE_B_SEQ, _PROBE_A_NAME, _PROBE_B_NAME, _TARGET_A_NAME, _TARGET_B_NAME,
            _TARGET_PAIR_ID, _PROBE_PAIR_ID]


def get_header_divider():
    return _HEADER_DIVIDER


def get_construct_header():
    return _CONSTRUCT_ID


def get_probe_id_header(probe_letter):
    return _PROBE_A_NAME if _is_letter_a(probe_letter) else _PROBE_B_NAME


def get_probe_seq_header(probe_letter):
    return _PROBE_A_SEQ if _is_letter_a(probe_letter) else _PROBE_B_SEQ


def get_target_id_header(target_letter):
    return _TARGET_A_NAME if _is_letter_a(target_letter) else _TARGET_B_NAME


def get_target_pair_id_header():
    return _TARGET_PAIR_ID


def get_probe_pair_id_header():
    return _PROBE_PAIR_ID


def get_comment_char():
    return "#"


def compose_probe_pair_id_from_probe_ids(probe_a_id, probe_b_id):
    divider = get_header_divider()
    result = probe_a_id + divider + divider + probe_b_id
    return result


def compose_target_pair_id_from_target_ids(target_a_id, target_b_id):
    return target_a_id + get_header_divider() + target_b_id


def extract_construct_and_grna_info(constructs_fp):
    construct_table = _read_in_construct_table(constructs_fp)
    seq_name_sets = _extract_unique_sets_across_a_and_b(construct_table,
        [_PROBE_A_NAME, _PROBE_A_SEQ], [_PROBE_B_NAME, _PROBE_B_SEQ])
    probe_name_seq_pairs = _validate_and_format_probe_seq_pairs(seq_name_sets)
    construct_names = construct_table[_CONSTRUCT_ID].unique().tolist()
    return construct_names, probe_name_seq_pairs


def trim_probes(probes_name_and_seq_list, retain_len):
    result = []
    for name_seq_tuple in probes_name_and_seq_list:
        probe_name = name_seq_tuple[0]
        full_seq = name_seq_tuple[1]
        trimmed_seq = trim_seq(full_seq, retain_len, False)  # False = do not retain from 5p end but from 3p end
        result.append((probe_name, trimmed_seq))
    return result


def load_annotation_df(constructs_fp, disregard_order):
    construct_df = _read_in_construct_table(constructs_fp)

    # TODO: I'm not thrilled that if there is non-uniqueness, I don't learn about it till this call;
    # I'd sort of like to learn about it from extract_construct_and_grna_info, even though that
    # method is able to successfully work around non-unique construct ids ...
    if len(construct_df[_CONSTRUCT_ID].unique()) != len(construct_df[_CONSTRUCT_ID]):
        raise ValueError("Non-unique construct ids detected.")

    if disregard_order:
        construct_df = _alphabetize_two_fields_in_row(construct_df,
                                                      get_target_id_header("a"), get_target_id_header("b"))
        construct_df = _alphabetize_two_fields_in_row(construct_df,
                                                      get_probe_id_header("a"), get_probe_id_header("b"))
    return construct_df


def _is_letter_a(letter):
    if letter.upper() == "A":
        result = True
    elif letter.upper() == "B":
        result = False
    else:
        raise ValueError("Input '{0}' is not recognized as A or B.".format(letter))

    return result


def _read_in_construct_table(constructs_fp, rows_to_skip=4):
    result = pandas.read_table(constructs_fp, comment=get_comment_char(), skiprows=rows_to_skip, header=None)
    result = _rename_columns(result)
    return result


def _rename_columns(construct_table, column_indices=None):
    new_names = [_CONSTRUCT_ID, _TARGET_A_NAME, _PROBE_A_NAME, _PROBE_A_SEQ,
                 _TARGET_B_NAME, _PROBE_B_NAME, _PROBE_B_SEQ]
    existing_names = list(construct_table.columns.values)
    if column_indices is None:
        column_indices = range(0, len(new_names))

    if len(column_indices) != len(new_names):
        raise ValueError("Expected indices for {0} columns but received indices for {1}.".format(
            len(new_names), len(column_indices)))

    existing_to_new_names = {}
    for curr_index in range(0, len(column_indices)):
        curr_col_index = column_indices[curr_index]
        curr_existing_name = existing_names[curr_col_index]
        existing_to_new_names[curr_existing_name] = new_names[curr_index]

    return construct_table.rename(columns=existing_to_new_names)


def _extract_unique_sets_across_a_and_b(construct_table, a_col_headers_list, b_col_headers_list):
    if len(a_col_headers_list) != len(b_col_headers_list):
        raise ValueError("A and B column header lists are not equal in length.")

    new_headers_list = ["temp_header_{0}".format(i) for i in range(0,len(a_col_headers_list))]

    # get the set of input columns for each of the two targets, assigning each the same
    # set of (generic) column headers so that they can easily be concatenated
    set_for_a = _extract_renamed_subset_df(construct_table, a_col_headers_list, new_headers_list)
    set_for_b = _extract_renamed_subset_df(construct_table, b_col_headers_list, new_headers_list)
    combined_set = pandas.concat([set_for_a, set_for_b])

    # extract only the unique sets
    grouped_combined_set = combined_set.groupby(new_headers_list).groups
    result = [x for x in grouped_combined_set]
    return sorted(result)  # NB sort so that output order is predictable


def _extract_renamed_subset_df(construct_table, col_headers_list, new_headers_list):
    result = construct_table.loc[:, col_headers_list]
    rename_dictionary = dict(zip(col_headers_list, new_headers_list))
    result.rename(columns=rename_dictionary, inplace=True)
    return result


def _validate_and_format_probe_seq_pairs(probes_seq_and_name_list):
    expected_num_pieces = 2
    seqs_by_names = {}
    names_by_seqs = {}
    result = []

    for curr_set in probes_seq_and_name_list:
        if len(curr_set) != expected_num_pieces:
            raise ValueError(
                "input '{0}' has {1} pieces instead of the expected {2}".format(
                    curr_set, len(curr_set), expected_num_pieces
                ))
        curr_seq = curr_set[1]
        curr_name = curr_set[0]

        if curr_seq in names_by_seqs:
            raise ValueError(
                "sequence '{0}' associated with name '{1}' but was already associated with name '{2}'".format(
                    curr_seq, curr_name, names_by_seqs[curr_seq]
                ))

        if curr_name in seqs_by_names:
            raise ValueError(
                "name '{0}' associated with sequence '{1}' but was already associated with sequence '{2}'".format(
                    curr_name, curr_seq, seqs_by_names[curr_name]
                ))

        names_by_seqs[curr_seq] = curr_name
        seqs_by_names[curr_name] = curr_seq

        result.append((curr_name, curr_seq.upper())) # upper-case all probe seqs
    # next pair in

    return result


def _alphabetize_two_fields_in_row(input_df, header_of_col_to_be_first, header_of_col_to_be_second):
    is_row_dealphabetized = input_df[header_of_col_to_be_first] >input_df[header_of_col_to_be_second] # boolean array

    # for dealphabetized rows, get header_of_col_to_be_first's value
    orig_first_col_vals_for_out_of_order_rows = input_df.loc[is_row_dealphabetized, header_of_col_to_be_first]

    # replace header_of_col_to_be_first's value with header_of_col_to_be_second's value for dealphabetized rows
    input_df.loc[is_row_dealphabetized, header_of_col_to_be_first] = input_df.loc[is_row_dealphabetized,
                                                                                header_of_col_to_be_second]
    # replace header_of_col_to_be_second's value with header_of_col_to_be_first's original value for dealphabetized rows
    input_df.loc[is_row_dealphabetized, header_of_col_to_be_second] = orig_first_col_vals_for_out_of_order_rows
    return input_df
