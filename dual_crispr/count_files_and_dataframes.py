# standard libraries
import os
import re

# third-party libraries
import pandas

# ccbb libraries
from ccbb_pyutils.analysis_run_prefixes import strip_run_prefix

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def get_counts_df(counts_fp, run_prefix):
    curr_counts_df = pandas.read_table(counts_fp, comment="#")  # TODO: remove hardcode
    orig_count_header = curr_counts_df.columns.values[-1]  # last column
    _, revised_count_header = os.path.split(orig_count_header)
    revised_count_header = strip_run_prefix(revised_count_header, run_prefix)
    curr_counts_df.rename(columns={orig_count_header: revised_count_header}, inplace=True)
    return revised_count_header, curr_counts_df


def clip_count_header_suffix(count_header):
    # if count header comes from IGM data out of Amanda's count pipeline, it will have
    # "_S#+_trimmed53_len_filtered_counts" on the end of it; get rid of this.
    # if count header comes out of Amanda's count pipeline, even if not from IGM, it will have
    # "_trimmed53_len_filtered_counts" on the end of it; get rid of this.
    # if it didn't come from these sources, it won't have this particular
    # suffix and the below trimming will simply have no effect.
    #
    # Regex key:
    # r means following is a raw string--many special characters ignored
    # [0-9]+ means "at least one digit, maybe more"
    # $ means "end of string"
    igm_pipeline_counts_regex = r'_S[0-9]+_trimmed53_len_filtered_counts$'
    pipeline_counts_regex = r'_trimmed53_len_filtered_counts$'
    result = re.sub(igm_pipeline_counts_regex, "", count_header)
    result = re.sub(pipeline_counts_regex, "", result)
    return result
