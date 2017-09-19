# third-party libraries
import pandas

# ccbb libraries
import dual_crispr.scoring_prep as ns_score_prep


class PerReplicateData(object):
    @staticmethod
    def get_mean_usable_log2_fraction_for_construct_header():
        return "mean_usable_log2_fraction_for_construct"

    @staticmethod
    def get_mean_of_usable_timepoints_header():
        return "mean_of_usable_timepoints"

    @staticmethod
    def get_variance_of_usable_timepoints_header():
        return "variance_of_usable_timepoints"

    @staticmethod
    def get_covariance_of_log2_fractions_w_timepoints_header():
        return "covariance_of_log2_fractions_w_timepoints"

    @staticmethod
    def extract_values_set_from_data_headers(data_only_df, get_timepts):
        values_set = set()
        data_headers = list(data_only_df.columns.values)
        for curr_data_header in data_headers:
            # read_timepoint_from_standardized_count_header--will already be cast to an integer
            curr_timept, curr_replicate = ns_score_prep.read_timepoint_and_replicate_from_standardized_count_header(
                curr_data_header)
            curr_value = curr_timept if get_timepts else curr_replicate
            values_set.add(curr_value)

        # get unique values and sort into ascending order
        return sorted(list(values_set))

    @property
    def replicate_num(self):
        return self._replicate_num

    @property
    def timepoints_series(self):
        return pandas.Series(data=self._timepoints_list,
                             index=self._log2_fractions_by_constructs_by_timepoints_df.columns.values)

    def __init__(self, replicate_num, log2_fractions_by_constructs_by_timepoints_df,
                 log2_fractions_thresholds_by_timepoints_series, time_prefixes_list=None):
        self._replicate_num = replicate_num
        self._log2_fractions_by_constructs_by_timepoints_df = log2_fractions_by_constructs_by_timepoints_df
        self._log2_fractions_thresholds_by_timepoints_series = log2_fractions_thresholds_by_timepoints_series
        if time_prefixes_list is None: time_prefixes_list = ["T","D"]
        self._time_prefixes_list = time_prefixes_list
        self._timepoints_list = self.extract_values_set_from_data_headers(
            self._log2_fractions_by_constructs_by_timepoints_df, get_timepts=True)
