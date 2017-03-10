# standard libraries
import re

# project-specific libraries
from ccbbucsd.malicrispr.construct_file_extracter import get_target_id_header, \
    get_probe_id_header, get_construct_header, get_target_pair_id_header, \
    get_probe_pair_id_header, compose_probe_pair_id_from_probe_ids, \
    compose_target_pair_id_from_target_ids, get_header_divider

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


_TIME_PREFIX = "T"
_NUM_HEADER_PIECES = 3


def get_time_prefix():
    return _TIME_PREFIX


def _get_num_header_pieces():
    return _NUM_HEADER_PIECES


def _clip_count_header_suffix(count_header):
    # if count header comes out of Amanda's count pipeline, it will have
    # "_S#+_trimmed53_len_filtered_counts" on the end of it; get rid of this.
    # if it didn't come out of Amanda's pipeline, it won't have this particular
    # suffix and the below trimming will simply have no effect.
    #
    # Regex key:
    # r means following is a raw string--many special characters ignored
    # [0-9]+ means "at least one digit, maybe more"
    # $ means "end of string"
    pipeline_counts_regex = r'_S[0-9]+_trimmed53_len_filtered_counts$'
    result = re.sub(pipeline_counts_regex, "", count_header)
    return result


def _validate_and_decompose_count_header(count_header):
    # Required count header format: experiment_timept_rep
    num_expected_pieces = _get_num_header_pieces()
    divider = get_header_divider()
    trimmed_count_header = _clip_count_header_suffix(count_header)

    count_header_pieces = trimmed_count_header.split(divider)
    return _validate_and_standardize_count_header_pieces(count_header_pieces)


def _validate_and_standardize_count_header_pieces(count_header_pieces):
    # Required count header format: experiment_timept_rep
    num_expected_pieces = _get_num_header_pieces()
    if len(count_header_pieces) != num_expected_pieces:
        raise ValueError("Count header has {0} pieces instead of the expected {1}: '{2}'.",
            len(count_header_pieces), num_expected_pieces, count_header_pieces)

    some_id = count_header_pieces[0]
    timept = _validate_and_standardize_timepoint(count_header_pieces[1])
    replicate = _validate_and_standardize_replicate(count_header_pieces[2])

    return (some_id, timept, replicate)


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
        raise ValueError("Time point value '{0}' is not recognizable as a positive integer.",
                         timept)

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


def validate_and_parse_data_column_headers(count_headers):
    # ensure that there is only one experiment represented in headers
    # ensure that every timepoint in the experiment has the exact same set of replicates

    expt_structure_by_id = {}
    result = []

    for curr_count_header in count_headers:
        # Required count header format: experiment_timept_rep
        count_header_pieces = _validate_and_decompose_count_header(curr_count_header)
        some_id = count_header_pieces[0]
        timept = count_header_pieces[1]
        replicate = count_header_pieces[2]
        result.append(count_header_pieces)

        # fill out structure {some_id: {timept: {set of replicates}}} for use in
        # validation after all columns are examined
        if some_id not in expt_structure_by_id: expt_structure_by_id[some_id] = {}
        curr_expt_structure = expt_structure_by_id[some_id]
        if timept not in curr_expt_structure: curr_expt_structure[timept] = set()
        curr_timept_replicates = curr_expt_structure[timept]
        curr_timept_replicates.add(replicate)

    _validate_expt_structure(expt_structure_by_id)
    return result


# So ... I wrote and tested _validate_expt_structure anticipating there could be
# MORE than ONE experiment in the input file.  Now I think that can't be the
# case.  However, it was a pain to write so in case I'm wrong, I'm leaving the
# version that handles that complexity here for now.
#
# def _validate_expt_structure(expt_structure_by_id):
#     # expt_structure_by_id should have format {some_id: {timept: {set of replicates}}}
#
#     # Every uid must have the exact same set of timepoints
#     # e.g., can't have sample1_T1_1, sample1_T2_1; sample2_T1_1, sample2_T2_1, sample2_T3_1
#     # Put all timepoints+replicates for a given sample in a set:
#     # set {T1_1, T2_1} doesn't match set {T1_1, T2_1, T3_1}
#
#     # All timepoints (across all expts) must have the exact same set of replicates.
#     # If timepts in *different* expts have different replicates, that would be
#     # caught by the above comparison of timept+replicate sets across samples.
#     # However, if timepts in the *same* sample have different numbers of replicates,
#     # that wouldn't be caught unless we keep sets by timept rather than timept+replicate.
#     # e.g., can't have sample1_T1_1; sample1_T2_1, sample1_T2_2
#
#     is_first_expt = True
#     reference_timepts_plus_reps_set = None
#
#     for curr_expt_id, curr_expt_structure in expt_structure_by_id.items():
#         # ensure all timepoints for current sample have the same number of replicates
#         is_first_timept = True
#         reference_reps_set = None
#
#         for curr_timept, curr_rep_set in curr_expt_structure.items():
#             if is_first_timept:
#                 reference_reps_set = curr_rep_set
#                 is_first_timept = False
#             else:
#                 if curr_rep_set != reference_reps_set:
#                     raise ValueError("For sample '{0}', timepoint {1} has "
#                                      "replicates '{2}' instead of the expected '{3}'".format(
#                         curr_expt_id, curr_timept, sorted(curr_rep_set),
#                         sorted(reference_reps_set)))
#
#             # make a new list of timept+rep for this timept
#
#
#         # TODO: handle case where no timepts exist, or no reps exist for timept
#         curr_timepts_plus_reps_set = {"{0}_{1}".format(timept, rep)
#                                  for timept in curr_expt_structure
#                                  for rep in reference_reps_set}
#
#         if is_first_expt:
#             reference_timepts_plus_reps_set = curr_timepts_plus_reps_set
#             is_first_expt = False
#         else:
#             if curr_timepts_plus_reps_set != reference_timepts_plus_reps_set:
#                 raise ValueError("Sample {0} has timepoints+replicates "
#                                  "'{1}' instead of the expected '{2}'".format(
#                     curr_expt_id, sorted(curr_timepts_plus_reps_set),
#                     sorted(reference_timepts_plus_reps_set)))
#
#
# def validate_data_column_headers(count_headers):
#     # ensure that every timepoint in an experiment has the exact same set of replicates
#     # ensure that every experiment has the exact same set of timepoints+replicates
#
#     expt_structure_by_id = {}
#     for curr_count_header in count_headers:
#         # Required count header format: experiment_timept_rep
#         count_header_pieces = _validate_and_decompose_count_header(curr_count_header)
#         some_id = count_header_pieces[0]
#         timept = _validate_and_standardize_timepoint(count_header_pieces[1])
#         replicate = count_header_pieces[2]
#
#         # fill out structure {some_id: {timept: {set of replicates}}} for use in
#         # validation after all columns are examined
#         if some_id not in expt_structure_by_id: expt_structure_by_id[some_id] = {}
#         curr_expt_structure = expt_structure_by_id[some_id]
#         if timept not in curr_expt_structure: curr_expt_structure[timept] = set()
#         curr_timept_replicates = curr_expt_structure[timept]
#         curr_timept_replicates.add(replicate)
#
#     _validate_expt_structure(expt_structure_by_id)


def _generate_scoring_friendly_annotation(annotation_df):
    construct_id_header = get_construct_header()
    target_pair_id_header = get_target_pair_id_header()
    probe_pair_id_header = get_probe_pair_id_header()

    divider = get_header_divider()

    result = annotation_df.loc[:, (construct_id_header, get_target_id_header("a"),
        get_probe_id_header("a"), get_target_id_header("b"), get_probe_id_header("b"))]
    target_pairs = (result[get_target_id_header("a")] + divider +
                    result[get_target_id_header("b")])
    result[target_pair_id_header] = result.apply(_compose_target_pair_id, axis=1)
    result[probe_pair_id_header] = result.apply(_compose_probe_pair_id, axis=1)
    return result


def _compose_probe_pair_id(row):
    return compose_probe_pair_id_from_probe_ids(row[get_probe_id_header("a")],
                                                row[get_probe_id_header("b")])


def _compose_target_pair_id(row):
    return compose_target_pair_id_from_target_ids(row[get_target_id_header("a")],
                                                row[get_target_id_header("b")])


def validate_and_recompose_count_header(expt_timeptnum_rep_tuple):
    standardized_pieces = _validate_and_standardize_count_header_pieces(expt_timeptnum_rep_tuple)

    divider = get_header_divider()
    result = "{}{}{}{}{}{}".format(standardized_pieces[0], divider, get_time_prefix(), standardized_pieces[1],
                divider, standardized_pieces[2])
    return result
