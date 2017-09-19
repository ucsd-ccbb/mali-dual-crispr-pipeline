# third-party libraries
import numpy

# ccbb libraries
import dual_crispr.construct_file_extracter as ns_extracter
import dual_crispr.scoring_prep as ns_score_prep
import dual_crispr.per_replicate_data as ns_per_rep

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class ScoringData:
    @staticmethod
    def _add_pseudocount(data_only_df, pseudocount):
        # NB: Added only to counts that are zero, not to all counts
        data_only_df[data_only_df == 0] = pseudocount
        return data_only_df

    @staticmethod
    def _calc_total_counts_per_sample(data_only_df):
        # axis = 0 means sum the values of the rows (over the column they are in)
        # columns are samples, so result is vector of 1 value per sample
        result = data_only_df.sum(axis=0)
        result.index = data_only_df.columns.values
        return result

    @staticmethod
    def _calc_log2_fractions(counts_for_row_construct_col_sample_df, total_counts_per_sample_series):
        # match the index of abundance_series to the columns of counts_for_row_construct_col_sample_df for dividing
        fraction_for_row_construct_col_sample_df = counts_for_row_construct_col_sample_df.div(
            total_counts_per_sample_series, axis="columns")
        log2_fraction_for_row_construct_col_sample_df = numpy.log2(fraction_for_row_construct_col_sample_df)
        return log2_fraction_for_row_construct_col_sample_df

    @staticmethod
    def _convert_abundance_thresholds_units(total_counts_per_sample_series, log2_counts_thresholds_per_sample_df):
        if not (log2_counts_thresholds_per_sample_df.index == total_counts_per_sample_series.index).all():
            raise ValueError(
                "Sample names for abundance thresholds do not match sample names for count data: [{0}], [{1}]".format(
                    ", ".join([str(x) for x in log2_counts_thresholds_per_sample_df.index]),
                    ", ".join([str(x) for x in total_counts_per_sample_series])
                ))

        log2_counts_thresholds_series = log2_counts_thresholds_per_sample_df.ix[:, 0]  # 0 = first column
        # subtraction in log space = division in "regular" space
        log2_fractions_thresholds_series = log2_counts_thresholds_series - numpy.log2(total_counts_per_sample_series)
        return log2_fractions_thresholds_series

    # NB: _run_init_methods param only exists to facilitate unit-testing and should NOT be changed for any other usage
    def __init__(self, scoring_prep_df, non_targeting_prefix, _run_init_methods=True):

        self._input_counts_and_annotation_df = scoring_prep_df
        self._non_targeting_prefix = non_targeting_prefix

        self._per_replicate_data_list = None
        self._genes = None
        self._probes = None

        if _run_init_methods:
            self._rename_non_targeting_gene()
            self._remove_constructs_w_both_probes_targeting_same_target_id()
            self._timepoints_list = ns_per_rep.PerReplicateData.extract_values_set_from_data_headers(
                self._get_data_only_df(), get_timepts=True)
            self._replicates_list = ns_per_rep.PerReplicateData.extract_values_set_from_data_headers(
                self._get_data_only_df(), get_timepts=False)

    @property
    def num_timepoints(self):
        return len(self._timepoints_list)

    @property
    def num_replicates(self):
        return len(self._replicates_list)

    # @property
    # def genes(self):
    #     if self._genes is None:
    #         self._genes = ns_extracter._extract_unique_sets_across_a_and_b(self._input_counts_and_annotation_df,
    #             [ns_extracter.get_target_id_header("A")], [ns_extracter.get_target_id_header("B")])
    #
    #     return self._genes
    #
    # @property
    # def probes(self):
    #     if self._probes is None:
    #         self._probes = ns_extracter._extract_unique_sets_across_a_and_b(self._input_counts_and_annotation_df,
    #             [ns_extracter.get_probe_id_header("A")], [ns_extracter.get_probe_id_header("B")])
    #
    #     return self._probes
    #
    # @property
    # def num_genes(self):
    #     return len(self.genes)

    def generate_per_replicate_data_list(self, log2_counts_thresholds_per_sample_df, pseudocount=1):
        data_only_df = self._get_data_only_df()
        pseudocounted_data_only_df = ScoringData._add_pseudocount(data_only_df, pseudocount)
        total_counts_per_sample_series = ScoringData._calc_total_counts_per_sample(pseudocounted_data_only_df)
        log2_fractions = ScoringData._calc_log2_fractions(pseudocounted_data_only_df, total_counts_per_sample_series)
        log2_fractions_thresholds_series = ScoringData._convert_abundance_thresholds_units(
            total_counts_per_sample_series, log2_counts_thresholds_per_sample_df)
        result = self._extract_replicate_data_list(log2_fractions, log2_fractions_thresholds_series)
        return result

    def _get_data_only_df(self):
        # NB: not guaranteed that all headers in this list are in dataframe
        annotation_headers = ns_extracter.get_potential_annotation_headers()
        non_annotation_headers = [x for x in self._input_counts_and_annotation_df if x not in annotation_headers]
        # NB: Must do a copy here so as to avoid potentially getting references to original df and changing accidentally
        data_only_df = self._input_counts_and_annotation_df[non_annotation_headers].copy()
        data_only_df.index = self._input_counts_and_annotation_df[ns_extracter.get_probe_pair_id_header()]
        return data_only_df

    def _remove_constructs_w_both_probes_targeting_same_target_id(self):
        not_same_target_id_mask = self._input_counts_and_annotation_df[ns_extracter.get_target_id_header("A")] != \
                                  self._input_counts_and_annotation_df[ns_extracter.get_target_id_header("B")]
        self._input_counts_and_annotation_df = self._input_counts_and_annotation_df[not_same_target_id_mask]

    def _rename_non_targeting_gene(self):
        # NB: It is acceptable to hard-code doing the renaming for both targets here, as there are only ever expected to
        # be TWO targets in a DUAL crispr screen.
        for curr_target_letter in ["A", "B"]:
            curr_target_col_header = ns_extracter.get_target_id_header(curr_target_letter)
            curr_target_series = self._input_counts_and_annotation_df[curr_target_col_header]
            non_targeting_target_rows_mask = curr_target_series.str.startswith(self._non_targeting_prefix)
            self._input_counts_and_annotation_df.loc[non_targeting_target_rows_mask, curr_target_col_header] = "0"

    def _get_column_headers_for_replicate(self, replicate_num):
        col_headers_for_replicate = []
        data_only_df = self._get_data_only_df()
        for curr_col_header in data_only_df.columns.values:
            _, curr_replicate_num = ns_score_prep.read_timepoint_and_replicate_from_standardized_count_header(
                curr_col_header)
            if curr_replicate_num == replicate_num:
                col_headers_for_replicate.append(curr_col_header)

        return col_headers_for_replicate

    def _extract_replicate_data_list(self, log2_fractions_df, log2_fractions_thresholds_series):
        replicate_data_list = []
        num_replicates = self.num_replicates
        for curr_replicate_num in range(1, num_replicates+1):  # replicate numbers start at 1, not 0
            relevant_sample_col_headers = self._get_column_headers_for_replicate(curr_replicate_num)
            log2_fractions_for_curr_replicate_df = log2_fractions_df[relevant_sample_col_headers].copy()
            log2_fractions_thresholds_for_curr_replicate_series = log2_fractions_thresholds_series[
                relevant_sample_col_headers].copy()
            curr_replicate_data = ns_per_rep.PerReplicateData(curr_replicate_num,
                                                              log2_fractions_for_curr_replicate_df,
                                                              log2_fractions_thresholds_for_curr_replicate_series)
            replicate_data_list.append(curr_replicate_data)

        return replicate_data_list



