# ccbb libraries
import ccbb_pyutils.pandas_utils as ns_pandas

# project-specific libraries
import dual_crispr.construct_file_extracter as ns_extracter
import dual_crispr.count_files_and_dataframes as ns_count

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


def get_prepped_file_suffix():
    return "timepoint_counts.txt"


def get_sample_name_header():
    return "sampleName"


def get_abundance_thresh_header():
    return "log2CountsThresh"


def get_abundance_thresh_file_suffix():
    return "abundance_thresholds.txt"


def read_timepoint_from_standardized_count_header(count_header, time_prefixes_list):
    count_header_pieces = count_header.split(ns_extracter.get_header_divider())
    if len(count_header_pieces) != _get_num_header_pieces() + 1:
        # +1 because we expect a standardized header, which has the experiment id added to it
        raise ValueError("Column header '{0}' splits into an unexpected number of pieces ({1})".format(
            count_header, len(count_header_pieces)
        ))

    timepoint_str = count_header_pieces[_get_timepoint_index()]
    return _validate_and_standardize_timepoint(timepoint_str, time_prefixes_list)


def merge_and_annotate_counts(count_file_fps, constructs_fp, dataset_name, time_prefixes_list,
                              disregard_order=True):
    construct_id_header = ns_extracter.get_construct_header()

    # load and merge the counts file(s)
    combined_counts_df = ns_pandas.merge_files_by_shared_header(count_file_fps, construct_id_header)

    # validate, standardize, and rename the count column headers
    orig_count_headers = _get_orig_count_headers(combined_counts_df)
    standardized_count_headers = _validate_and_standardize_count_headers(orig_count_headers,
                                                                         dataset_name, time_prefixes_list)
    rename_dictionary = dict(zip(orig_count_headers, standardized_count_headers))
    combined_counts_df.rename(columns=rename_dictionary, inplace=True)

    # load and standardize the annotation file (containing construct definitions)
    annotation_df = ns_extracter.load_annotation_df(constructs_fp, disregard_order)
    minimal_annotation_df = _generate_scoring_friendly_annotation(annotation_df)

    # join counts to annotation and sort into required order
    joined_df = minimal_annotation_df.merge(combined_counts_df, on=construct_id_header)
    sorted_col_headers = list(minimal_annotation_df.columns.values)
    sorted_count_headers = _sort_headers(standardized_count_headers, dataset_name, time_prefixes_list)
    sorted_col_headers.extend(sorted_count_headers)
    return joined_df.loc[:, sorted_col_headers]


def _get_timepoint_index():
    return -2


def _get_num_header_pieces():
    return 2


def _get_preferred_timept_prefix(time_prefixes_list):
    return time_prefixes_list[0]


def _get_orig_count_headers(combined_counts_df):
    orig_count_headers = list(combined_counts_df.columns.values)
    orig_count_headers.remove(ns_extracter.get_construct_header())
    return orig_count_headers


def _validate_and_standardize_count_headers(orig_count_headers, expt_id, time_prefixes_list):
    error_msgs = []
    expt_structure_by_id = {}
    orig_header_by_standardized = {}
    result = []

    for curr_count_header in orig_count_headers:
        # Required count header format: experiment_timept_rep
        try:
            valid_id, timept, replicate = _validate_and_standardize_header_pieces(curr_count_header, expt_id,
                                                                                  time_prefixes_list)

            standardized_header = _recompose_count_header(valid_id, timept, replicate, time_prefixes_list)
            if standardized_header in orig_header_by_standardized:
                raise ValueError("The following pair of column headers both appear to represent the same timepoint "
                                 "and replicate: '{0}', '{1}'.  Please modify the inputs to remove this "
                                 "ambiguity.".format(curr_count_header, orig_header_by_standardized[
                                 standardized_header]))
        except ValueError as ex:
            error_msgs.append(str(ex))
            continue

        orig_header_by_standardized[standardized_header] = curr_count_header
        result.append(standardized_header)

        # fill out structure {some_id: {timept: {set of replicates}}} for use in
        # validation after all columns are examined
        if valid_id not in expt_structure_by_id: expt_structure_by_id[valid_id] = {}
        curr_expt_structure = expt_structure_by_id[valid_id]
        if timept not in curr_expt_structure: curr_expt_structure[timept] = set()
        curr_timept_replicates = curr_expt_structure[timept]
        curr_timept_replicates.add(replicate)

    if len(error_msgs) > 0:
        raise ValueError(
            'The following error(s) were detected during count file parsing:\n{0}'.format("\n".join(error_msgs)))

    _validate_expt_structure(expt_structure_by_id)
    return result


# this method is broken out from _validate_and_standardize_count_headers just to make unit testing easier
def _validate_and_standardize_header_pieces(curr_count_header, expt_id, time_prefixes_list):
    timept, replicate = _validate_and_standardize_timept_and_replicate(
        curr_count_header, time_prefixes_list)
    valid_id = _validate_and_standardize_expt_id(expt_id)

    return valid_id, timept, replicate


def _validate_and_standardize_timept_and_replicate(count_header, time_prefixes_list):
    # Required count header format: experiment_timept_rep
    trimmed_count_header = ns_count.clip_count_header_suffix(count_header)
    count_header_pieces = trimmed_count_header.split(ns_extracter.get_header_divider())

    num_expected_pieces = _get_num_header_pieces()
    if len(count_header_pieces) < num_expected_pieces:
        raise ValueError("Column header '{0}' separates on the '{1}' delimiter into the following {2} piece(s) instead of the expected {3}: "
                         "{3}.".format(count_header, ns_extracter.get_header_divider(), len(count_header_pieces),
                                       num_expected_pieces, count_header_pieces))

    timept = _validate_and_standardize_timepoint(count_header_pieces[_get_timepoint_index()],
                                                 time_prefixes_list)
    replicate = _validate_and_standardize_replicate(count_header_pieces[-1])
    return timept, replicate


def _validate_and_standardize_timepoint(timept, time_prefixes_list):
    if isinstance(timept, str):
        # ensure timepoint is "t" or "T" plus a non-negative integer number
        expected_timepoint_prefixes = [x.upper() for x in time_prefixes_list]
        timepoint_prefix = timept[:1]
        if timepoint_prefix.upper() not in expected_timepoint_prefixes:
            raise ValueError("Time point '{0}' does not start with upper or lower case versions of any of the "
                             "expected prefixes {1}.".format(
                timept, ", ".join(time_prefixes_list)))

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


def _recompose_count_header(expt_id, timept, replicate, time_prefixes_list):
    time_prefix = _get_preferred_timept_prefix(time_prefixes_list)
    divider = ns_extracter.get_header_divider()
    result = "{}{}{}{}{}{}".format(expt_id, divider, time_prefix,
                                   timept, divider, replicate)
    return result


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
    # target_pair_id_header = ns_extracter.get_target_pair_id_header()
    # probe_pair_id_header = ns_extracter.get_probe_pair_id_header()
    # Note: the below column creations could be done without using apply (i.e., by
    # just writing "df[colA] + divider + df[colB]") but I used apply because I want
    # to centralize the code that creates these strings, and sometimes it needs to work
    # on a single pair of variables rather than columns of variables, so it needed to be
    # a non-vectorized function.
    # result[target_pair_id_header] = result.apply(_compose_target_pair_id, axis=1)
    # result[probe_pair_id_header] = result.apply(_compose_probe_pair_id, axis=1)
    return result

# def _compose_probe_pair_id(row):
#     return ns_extracter.compose_probe_pair_id_from_probe_ids(row[ns_extracter.get_probe_id_header("a")],
#                                                 row[ns_extracter.get_probe_id_header("b")])
#
#
# def _compose_target_pair_id(row):
#     return ns_extracter.compose_target_pair_id_from_target_ids(row[ns_extracter.get_target_id_header("a")],
#                                                   row[ns_extracter.get_target_id_header("b")])


def _sort_headers(headers_list, expt_id, time_prefixes_list):
    # headers need to be sorted first by timepoint *as a number* and then by replicate *as a number*
    # the expt id should be irrelevant to the sorting as it should be the same for all headers
    header_tuples_list = [_validate_and_standardize_header_pieces(x, expt_id,time_prefixes_list) for x in headers_list]
    sorted_header_tuples_list = sorted(header_tuples_list)
    result = [_recompose_count_header(*x, time_prefixes_list=time_prefixes_list) for x in sorted_header_tuples_list]
    return result
