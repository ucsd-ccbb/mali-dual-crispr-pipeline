# third-party libraries
import numpy
import pandas

# ccbb libraries
import dual_crispr.construct_file_extracter as ns_extracter
import dual_crispr.scoring_prep as ns_score_prep
import dual_crispr.per_replicate_data as ns_per_rep

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "development"


class ScoringData:
    # Method generalizes functionality from Roman's MethodII.R line 391, shown below with additional comments by ab:
    # ab: gooddata = dataframe of input counts per construct (row) by sample (column), scrubbed of any rows where both
    # ab: probes in construct target same gene.
    # gooddata[gooddata == 0] <- 1  # pseudocounts
    # Note that unlike R code, this code parameterizes pseudocount so it can be changed if desired.
    @staticmethod
    def _add_pseudocount(data_only_df, pseudocount):
        # NB: Added only to counts that are zero, not to all counts
        data_only_df[data_only_df == 0] = pseudocount
        return data_only_df  # i.e., gooddata after pseudocounting

    # Method implements functionality from Roman's MethodII.R line 392:
    # ab: gooddata = dataframe of input counts per construct (row) by sample (column), scrubbed of any rows where both
    # ab: probes in construct target same gene and with any 0 cells given a pseudocount of 1.
    # abundance<-apply(gooddata,2,sum)
    @staticmethod
    def _calc_total_counts_per_sample(data_only_df):  # i.e., gooddata after pseudocounting
        # axis = 0 means sum the values of the rows (over the column they are in)
        # columns are samples, so result is vector of 1 value per sample
        result = data_only_df.sum(axis=0)
        result.index = data_only_df.columns.values
        return result  # i.e., abundance

    # Method extends Roman's MethodII.R line 393:
    # ab: gooddata = dataframe of input counts per construct (row) by sample (column), scrubbed of any rows where both
    # ab: probes in construct target same gene and with any 0 cells given a pseudocount of 1.
    # ab: abundance = total number of counts in gooddata per sample
    # y<-t(log2(t(gooddata)/abundance)) #log2 frequencies
    @staticmethod
    def _calc_log2_fractions(counts_for_row_construct_col_sample_df, total_counts_per_sample_series):
        # i.e., inputs are gooddata and abundance.
        # match the index of abundance_series to the columns of counts_for_row_construct_col_sample_df for dividing
        fraction_for_row_construct_col_sample_df = counts_for_row_construct_col_sample_df.div(
            total_counts_per_sample_series, axis="columns")
        log2_fraction_for_row_construct_col_sample_df = numpy.log2(fraction_for_row_construct_col_sample_df)
        return log2_fraction_for_row_construct_col_sample_df  # i.e., y

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

        if _run_init_methods:
            # TODO: Think about whether I need to use the "rename non-targeting gene to 0" strategy or
            # whether I can come up with something less hacky
            # self._rename_non_targeting_gene()

            self._remove_constructs_w_both_probes_targeting_same_target_id()
            self._input_counts_and_annotation_df.set_index([ns_extracter.get_probe_pair_id_header()], inplace=True)
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

    # TODO: Decide if I can get rid of these ... not used?
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

    def annotate_per_construct_series_to_df(self, per_construct_series):
        gene_probe_annotation_only_df = self._input_counts_and_annotation_df[[ns_extracter.get_target_id_header("A"),
                                                                              ns_extracter.get_probe_id_header("A"),
                                                                              ns_extracter.get_target_id_header("B"),
                                                                              ns_extracter.get_probe_id_header("B")]]
        annotated_per_construct_df = pandas.concat([per_construct_series.to_series, gene_probe_annotation_only_df],
                                                   axis=1)
        return annotated_per_construct_df

    def make_by_tuple_by_tuple_annotated_symmetric_df(self, per_construct_series):
        annotated_per_construct_df = self.annotate_per_construct_series_to_df(per_construct_series)
        geneA_header = ns_extracter.get_target_id_header("A")
        geneB_header = ns_extracter.get_target_id_header("B")
        probeA_header = ns_extracter.get_probe_id_header("A")
        probeB_header = ns_extracter.get_probe_id_header("B")
        gene_probe_A_tuple_header = ns_extracter.get_gene_probe_tuple_header("A")
        gene_probe_B_tuple_header = ns_extracter.get_gene_probe_tuple_header("B")

        # add new columns containing tuple of gene, probe for target A and target B of construct
        annotated_per_construct_df[gene_probe_A_tuple_header] = list(zip(annotated_per_construct_df[geneA_header],
                                                   annotated_per_construct_df[probeA_header]))
        annotated_per_construct_df[gene_probe_B_tuple_header] = list(zip(annotated_per_construct_df[geneB_header],
                                                   annotated_per_construct_df[probeB_header]))

        # now remove the columns that I made the tuples from
        annotated_per_construct_df = annotated_per_construct_df.drop([geneA_header, probeA_header, geneB_header,
                                                                      probeB_header], axis=1)

        # data values are in the first column now, as all columns that came before that in df have been dropped
        values_header = annotated_per_construct_df.columns.values[0]

        # pivot the dataframe so it becomes matrix-like, with tuples for target a as row labels, tuples for target b
        # as column labels, and data values in each cell. Note it isn't yet symmetric--may have *different* tuple
        # values on each axis.
        annotated_by_tuple_by_tuple_df = annotated_per_construct_df.pivot(index=gene_probe_A_tuple_header,
                                                                          columns=gene_probe_B_tuple_header,
                                                                          values=values_header)

        # now transpose the matrix-like df and add the transpose to itself.  This creates a symmetric matrix-like df
        # with all potential values in both sets of tuples occurring (in the same order) on both axes.  In cases where
        # the addition would be of an NaN plus a number, pretend the NaN is 0 (hence, fill_value=0).  Note that when
        # BOTH cells being added are NaN, the output will STILL be NaN rather than 0+0=0, so see next line.
        symmetric_by_tuple_by_tuple_df = annotated_by_tuple_by_tuple_df.add(annotated_by_tuple_by_tuple_df.transpose(),
                                                                            fill_value=0)
        # Since some tuples have no data (e.g., we never pair a probe with itself in a construct, or with another probe
        # for its own gene), fill in those NaNs with 0.  Why 0? Because Roman's original code uses 0 as the base fill
        # for all its symmetric matrices.
        symmetric_by_tuple_by_tuple_df.fillna(0, inplace=True)
        return symmetric_by_tuple_by_tuple_df

    # Method extends Roman's MethodII.R lines 390-393:
    # ab: gooddata = dataframe of annotations for each construct and input counts per construct (row) by sample
    # ab: (column), scrubbed of any rows where both probes in construct target same gene.
    #
    # gooddata<-data.matrix(goodX[,6:(5+2*nt)])
    # gooddata[gooddata==0]<-1 #pseudocounts
    # abundance<-apply(gooddata,2,sum)
    # y<-t(log2(t(gooddata)/abundance)) #log2 frequencies
    #
    # Method also produces equivalent of "ab0" variable from Roman's MethodII.R line 5:
    # ab: ab0 = log2 of counts thresholds by sample
    # ab0<-c(-18.,-17.5,-18.,-18.,-18,-18,-17.5,-17.5)
    # from input data that are in different units.
    # Note that unlike R code, this code doesn't assume number of annotation columns or depend on number of timepoints.
    def generate_per_replicate_data_list(self, log2_counts_thresholds_per_sample_df, pseudocount=1):
        data_only_df = self._get_data_only_df()  # i.e., gooddata after slicing
        pseudocounted_data_only_df = ScoringData._add_pseudocount(data_only_df, pseudocount)
        # total_counts_per_sample_series = abundance
        total_counts_per_sample_series = ScoringData._calc_total_counts_per_sample(pseudocounted_data_only_df)
        # log2_fractions = y
        log2_fractions = ScoringData._calc_log2_fractions(pseudocounted_data_only_df, total_counts_per_sample_series)
        # log2_fractions_thresholds_series = ab0
        log2_fractions_thresholds_series = ScoringData._convert_abundance_thresholds_units(
            total_counts_per_sample_series, log2_counts_thresholds_per_sample_df)
        result = self._extract_replicate_data_list(log2_fractions, log2_fractions_thresholds_series)
        return result

    # Method generalizes dataframe slice X[,6:(5+2*nt)] from Roman's MethodII.R (e.g., line 354, 390), where was
    # expected that first 5 columns are annotation, data starts in 6; nt = num timepts; 2 = number of replicates.
    # Note that unlike R code, this code doesn't assume number of annotation columns or depend on number of timepoints.
    def _get_data_only_df(self):
        # NB: not guaranteed that all headers in this list are in dataframe
        annotation_headers = ns_extracter.get_potential_annotation_headers()
        non_annotation_headers = [x for x in self._input_counts_and_annotation_df if x not in annotation_headers]
        # NB: Must do a copy here so as to avoid potentially getting references to original df and changing accidentally
        data_only_df = self._input_counts_and_annotation_df[non_annotation_headers].copy()
        data_only_df.index = self._input_counts_and_annotation_df[ns_extracter.get_probe_pair_id_header()]
        return data_only_df

    # Method implements Roman's MethodII.R line 355:
    # ab: X<-read.table(input_filename,sep="\t",header=TRUE)
    # ab: geneA and geneB are hard-coded expected column names
    # good <- (X$geneA != X$geneB)
    # Note that this method DOES still assume fixed column names, although those aren't hard-coded here.
    def _remove_constructs_w_both_probes_targeting_same_target_id(self):
        not_same_target_id_mask = self._input_counts_and_annotation_df[ns_extracter.get_target_id_header("A")] != \
                                  self._input_counts_and_annotation_df[ns_extracter.get_target_id_header("B")]
        self._input_counts_and_annotation_df = self._input_counts_and_annotation_df[not_same_target_id_mask]

    # def _rename_non_targeting_gene(self):
    #     # NB: It is acceptable to hard-code doing the renaming for both targets here, as only ever expect
    #     # TWO targets in a DUAL crispr screen.
    #     for curr_target_letter in ["A", "B"]:
    #         curr_target_col_header = ns_extracter.get_target_id_header(curr_target_letter)
    #         curr_target_series = self._input_counts_and_annotation_df[curr_target_col_header]
    #         non_targeting_target_rows_mask = curr_target_series.str.startswith(self._non_targeting_prefix)
    #         self._input_counts_and_annotation_df.loc[non_targeting_target_rows_mask, curr_target_col_header] = "0"

    def _get_column_headers_for_replicate(self, replicate_num):
        col_headers_for_replicate = []
        data_only_df = self._get_data_only_df()
        for curr_col_header in data_only_df.columns.values:
            _, curr_replicate_num = ns_score_prep.read_timepoint_and_replicate_from_standardized_count_header(
                curr_col_header)
            if curr_replicate_num == replicate_num:
                col_headers_for_replicate.append(curr_col_header)

        return col_headers_for_replicate

    # Method extends Roman's MethodII.R lines 418-419:
    # ab: y = log2_fractions per construct per sample
    # ab: ab0 = log2 of counts thresholds by sample
    # ab: nt = number of timepoints
    # ab: 2 is assumed number of replicates
    # ab: assumption is that sample columns are ordered by timepoint and within timepoint by replicate--
    # ab: e.g. myexpt_timept1_rep1  myexpt_timept1_rep2 myexpt_timept2_rep1 myexpt_timept2_rep2--
    # ab: so code is pulling out columns samples for all timepoints for a single replicate
    # x1 <- y[, seq(1, 2 * nt, by=2)]
    # ab1 <- ab0[seq(1, 2 * nt, by=2)]
    # Note that unlike R code, this method doesn't assume two replicates and doesn't assume order of columns.
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



