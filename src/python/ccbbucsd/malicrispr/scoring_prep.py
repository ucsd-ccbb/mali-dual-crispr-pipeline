# ccbb libraries
import ccbbucsd.utilities.pandas_utils as ns_pandas

# project-specific libraries
import ccbbucsd.malicrispr.construct_file_extracter as ns_extracter
import ccbbucsd.malicrispr.count_files_and_dataframes as ns_count


__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"

_NUM_HEADER_PIECES = 3


def get_time_prefix():
    return "T"


def get_prepped_file_suffix():
    return "timepoint_counts.txt"


def get_sample_name_header():
    return "sampleName"


def get_abundance_thresh_header():
    return "log2CountsThresh"


def get_abundance_thresh_file_suffix():
    return "abundance_thresholds.txt"


def merge_and_annotate_counts(count_file_fps, constructs_fp, column_indices,
                              dataset_name, disregard_order=True):

    construct_id_header = ns_extracter.get_construct_header()

    # load and validate the counts file(s)
    combined_counts_df = ns_pandas.merge_files_by_shared_header(count_file_fps, construct_id_header)
    header_pieces_tuples_list = _validate_and_rename_counts_columns(combined_counts_df, dataset_name)
    #expt_name = header_pieces_tuples_list[0][0]

    # load and standardize the annotation file (containing construct definitions)
    annotation_df = ns_extracter.load_annotation_df(constructs_fp, column_indices, disregard_order)
    minimal_annotation_df = _generate_scoring_friendly_annotation(annotation_df)

    # join counts to annotation and sort into required order
    joined_df = minimal_annotation_df.merge(combined_counts_df, on=construct_id_header)
    sorted_col_headers = _get_sorted_joined_df_column_headers(minimal_annotation_df,
                                                              header_pieces_tuples_list)
    return dataset_name, joined_df.loc[:, sorted_col_headers]


def _validate_and_rename_counts_columns(combined_counts_df, expt_id):
    data_col_headers = list(combined_counts_df.columns.values)
    data_col_headers.remove(ns_extracter.get_construct_header())

    result = _validate_and_parse_data_column_headers(data_col_headers, expt_id)

    # recompose headers rather than using data_col_headers because their components
    # may have been changed by standardization done in _validate_and_parse_data_column_headers
    unsorted_headers_list = _recompose_headers_from_tuples(result)
    rename_dictionary = dict(zip(data_col_headers, unsorted_headers_list))
    combined_counts_df.rename(columns=rename_dictionary, inplace=True)

    return result


def _validate_and_parse_data_column_headers(count_headers, expt_id):
    # ensure that there is only one experiment represented in headers
    # ensure that every timepoint in the experiment has the exact same set of replicates

    expt_structure_by_id = {}
    result = []

    for curr_count_header in count_headers:
        # Required count header format: experiment_timept_rep
        count_header_pieces = _validate_and_decompose_count_header(curr_count_header)
        some_id = _validate_and_standardize_expt_id(expt_id) # count_header_pieces[0]
        timept = count_header_pieces[1]
        replicate = count_header_pieces[2]
        result.append((some_id, timept, replicate))  # TODO: Ugly hacking here; come back and clean up

        # fill out structure {some_id: {timept: {set of replicates}}} for use in
        # validation after all columns are examined
        if some_id not in expt_structure_by_id: expt_structure_by_id[some_id] = {}
        curr_expt_structure = expt_structure_by_id[some_id]
        if timept not in curr_expt_structure: curr_expt_structure[timept] = set()
        curr_timept_replicates = curr_expt_structure[timept]
        curr_timept_replicates.add(replicate)

    _validate_expt_structure(expt_structure_by_id)
    return result


def _validate_and_decompose_count_header(count_header):
    # Required count header format: experiment_timept_rep
    trimmed_count_header = ns_count.clip_count_header_suffix(count_header)
    count_header_pieces = trimmed_count_header.split(ns_extracter.get_header_divider())
    return _validate_and_standardize_count_header_pieces(count_header_pieces)


def _validate_and_standardize_count_header_pieces(count_header_pieces):
    num_expected_pieces = _get_num_header_pieces()
    if len(count_header_pieces) != num_expected_pieces:
        raise ValueError("Count header has {0} piece(s) instead of the expected {1}: '{2}'.".format(
                         len(count_header_pieces), num_expected_pieces, count_header_pieces))

    some_id = count_header_pieces[0]  # _validate_and_standardize_expt_id(count_header_pieces[0])
    timept = _validate_and_standardize_timepoint(count_header_pieces[1])
    replicate = _validate_and_standardize_replicate(count_header_pieces[2])

    return some_id, timept, replicate


def _get_num_header_pieces():
    return _NUM_HEADER_PIECES


def _validate_and_standardize_expt_id(expt_id):
    # experiment ids are used as a component of the sample names--and sample names are used as
    # dataframe *column names* in both pandas (Python) and R.  In R, column names must be valid (R) variable
    # names, and thus may contain only alphanumerics, periods, and underscores.  In pandas, column names
    # must be valid Python variable names, which may contain only alphanumerics and underscores.  Thus,
    # the only *non*-alphanumeric character that is accepted by both is the underscore, and I use that
    # to delimit the pieces of the sample name (i.e., exptid_timept_replicatenum).  That means expt id
    # can't contain any underscores itself--leaving only alphanumerics :(
    if not expt_id.isalnum():
        raise ValueError("Experiment id '{0}' is not strictly alphanumeric.".format(
                         expt_id))

    return expt_id


def _validate_and_standardize_timepoint(timept):
    if isinstance(timept, str):
        # ensure timepoint is "t" or "T" plus a non-negative integer number
        expected_timepoint_prefix = get_time_prefix()
        timepoint_prefix = timept[:1]
        if timepoint_prefix.upper() != expected_timepoint_prefix.upper():
            raise ValueError("Time point '{0}' does not start with '{1}' or '{2}'.", timept,
                             expected_timepoint_prefix.lower(), expected_timepoint_prefix.upper())

        timept = timept[1:]
    else:
        timept = str(timept)

    if not timept.isdigit():
        raise ValueError("Time point value '{0}' is not recognizable as a positive integer.".format(
                         timept))

    return int(timept)


def _validate_and_standardize_replicate(rep):
    if not isinstance(rep, int):
        rep = int(rep) if rep.isdigit() else rep
    return rep


def _validate_expt_structure(expt_structure_by_id):
    # expt_structure_by_id should have format {some_id: {timept: {set of replicates}}}

    # There must be only one experiment represented in the data structure

    # All timepoints in the experiment must have the exact same set of replicates:
    # e.g., can't have sample1_T1_1; sample1_T2_1, sample1_T2_2

    if len(expt_structure_by_id) != 1:
        raise ValueError(("Count headers must describe one and only one experiment, "
                          "but {0} were detected: '{1}'.").format(len(expt_structure_by_id),
                                                                  sorted(list(expt_structure_by_id.keys()))))

    for curr_expt_id, curr_expt_structure in expt_structure_by_id.items():
        # ensure all timepoints for current sample have the same number of replicates
        is_first_timept = True
        reference_reps_set = None

        if len(curr_expt_structure) == 0:
            raise ValueError("Count headers must describe at least one timepoint for experiment, "
                             "but 0 were detected.")

        for curr_timept, curr_rep_set in curr_expt_structure.items():
            if len(curr_rep_set) == 0:
                raise ValueError(("Count headers must describe at least one replicate for each timepoint, "
                                  "but 0 were detected for timepoint '{0}'.").format(curr_timept))

            if is_first_timept:
                reference_reps_set = curr_rep_set
                is_first_timept = False
            else:
                if curr_rep_set != reference_reps_set:
                    raise ValueError("For sample '{0}', timepoint {1} has "
                                     "replicates '{2}' instead of the expected '{3}'".format(
                        curr_expt_id, curr_timept, sorted(curr_rep_set),
                        sorted(reference_reps_set)))


def _recompose_headers_from_tuples(header_pieces_tuples_list, sort=False):
    input_list = sorted(header_pieces_tuples_list) if sort else header_pieces_tuples_list
    return [_validate_and_recompose_count_header(x) for x in input_list]


def _generate_scoring_friendly_annotation(annotation_df):
    construct_id_header = ns_extracter.get_construct_header()

    result = annotation_df.loc[:, [construct_id_header,
                                   ns_extracter.get_probe_id_header("a"),
                                   ns_extracter.get_probe_id_header("b"),
                                   ns_extracter.get_target_id_header("a"),
                                   ns_extracter.get_target_id_header("b")
                                   ]]

    # Below is what I expect the output to be after scoring data prep code is refactored to accept more
    # detail (and generate less itself).
    # result = annotation_df.loc[:, (construct_id_header, ns_extracter.get_target_id_header("a"),
    #                                ns_extracter.get_probe_id_header("a"),
    #                                ns_extracter.get_target_id_header("b"),
    #                                ns_extracter.get_probe_id_header("b"))]
    #target_pair_id_header = ns_extracter.get_target_pair_id_header()
    #probe_pair_id_header = ns_extracter.get_probe_pair_id_header()
    # Note: the below column creations could be done without using apply (i.e., by
    # just writing "df[colA] + divider + df[colB]") but I used apply because I want
    # to centralize the code that creates these strings, and sometimes it needs to work
    # on a single pair of variables rather than columns of variables, so it needed to be
    # a non-vectorized function.
    #result[target_pair_id_header] = result.apply(_compose_target_pair_id, axis=1)
    #result[probe_pair_id_header] = result.apply(_compose_probe_pair_id, axis=1)
    return result


# def _compose_probe_pair_id(row):
#     return ns_extracter.compose_probe_pair_id_from_probe_ids(row[ns_extracter.get_probe_id_header("a")],
#                                                 row[ns_extracter.get_probe_id_header("b")])
#
#
# def _compose_target_pair_id(row):
#     return ns_extracter.compose_target_pair_id_from_target_ids(row[ns_extracter.get_target_id_header("a")],
#                                                   row[ns_extracter.get_target_id_header("b")])


def _validate_and_recompose_count_header(expt_timeptnum_rep_tuple):
    standardized_pieces = _validate_and_standardize_count_header_pieces(expt_timeptnum_rep_tuple)

    divider = ns_extracter.get_header_divider()
    result = "{}{}{}{}{}{}".format(standardized_pieces[0], divider, get_time_prefix(), standardized_pieces[1],
                                   divider, standardized_pieces[2])
    return result


def _get_sorted_joined_df_column_headers(minimal_annotation_df, header_pieces_tuples_list):
    result = list(minimal_annotation_df.columns.values)
    sorted_data_headers = _recompose_headers_from_tuples(header_pieces_tuples_list, True)
    result.extend(sorted_data_headers)
    return result
