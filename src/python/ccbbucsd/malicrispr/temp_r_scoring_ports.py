# third-party libraries
import pandas

# return a data frame containing only rows whose values in
# the given column are in the input list, but all columns.
# if a column name is not specified, return only rows whose
# row names are in the input list.
def get_subset_of_df_rows(input_df, input_values, col_name_to_check=None, exclude_input_values=False):
    row_values_to_check = list(input_df.index) if col_name_to_check is None else input_df[[col_name_to_check]]
    rows_to_include = [x for x in row_values_to_check if x in input_values]
    if exclude_input_values:
        rows_to_include = not rows_to_include

    result = input_df[rows_to_include,]
    return result


# class ReplicateDataStore():
#     def __init__(self, full_screen_df, replicate_num):
#         self._construct_annotation = None
#         self._log2freq_matrix = None
#         self._passes_thresh_matrix = None
#
#     def _generate_threshold_passing_matrix(self, abundance_threshs_by_sample_df):
#
#
#
#
# # log2FreqsConstructBySampleMatrix should have constructs as rows and samples for this replicate as columns.
# # abundanceThresholdsBySampleDf is expected to have construct ids as rownames.
# generateThresholdPassingMatrix < - function(log2FreqsConstructBySampleMatrix, abundanceThresholdsBySampleDf)
# {
#     # get names of all samples for this replicate from the column names of aReplicateData$log2FreqsMatrix
#     # then subset abundanceThresholdBySampleDf by removing any rows whose rowname is not in that sample list
#     relevantAbundanceThresholdsDf = getSubsetOfDfRows(abundanceThresholdsBySampleDf,
#                                                       colnames(log2FreqsConstructBySampleMatrix),
#                                                       colNameToCheck=NULL,
#                                                       excludeInputValues=FALSE)
#
# # transpose the log2FreqsMatrixForReplicate so samples are rows and constructs are in columns, and cell values are log2frequencies
# # generate new matrix containing w/same rows and columns, but cell values are whether that construct+sample passed the sample's abundance threshold
# # transpose new matrix so that constructs are in rows and samples are in columns
# log2FreqsSampleByConstructMatrix = t(log2FreqsConstructBySampleMatrix)
# passesAbundanceThreshSampleByConstructMatrix = log2FreqsSampleByConstructMatrix > relevantAbundanceThresholdsDf
# passesAbundanceThreshConstructBySampleMatrix = t(passesAbundanceThreshSampleByConstructMatrix)
# return (passesAbundanceThreshConstructBySampleMatrix)
# }
#
#
#
#
#
#
#
# # constructs as rows, timept for 1st replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2 freq for row/col combination is above relevant abundance threshold
# good1 < - t(t(x1) > ab1)
# # constructs as rows, timept for 2nd replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2 freq for row/col combination is above relevant abundance threshold
# good2 < - t(t(x2) > ab2)
# # 1 = sum over rows--i.e., constructs.  Any construct that isn't above abundance threshold in 1st replicate in at least 2 timepoints has 1/TRUE in "useless" matrix
# useless1 < - apply(good1, 1, sum) < 2
# useless2 < - apply(good2, 1, sum) < 2
#
# # ab
# # print(sum(useless1))
# # print(sum(useless2))
# # PREDICTION CORRECT!
# # --------------------
#
# # I do not understand how the statements below actually *remove* anything from the good1 and good2 dfs ...
# # It seems like, for each constructs that isn't above abundance thresholds in at least two timepoints for this replicate, it sets that construct's "good" values to FALSE for *all* timepoints in this replicate
# good1[useless1,] < - FALSE  # remove singletons
# good2[useless2,] < - FALSE  # remove singletons
#
# }
#
#
# idConstructsFailingAbudanceThresholds < - function(aReplicateObj, abundanceThresholdsBySampleDf,
#                                                    minNumSamplesAboveAbundance)
# {
#     passesAbundanceThreshConstructBySampleMatrix = generateThresholdPassingMatrix(
#     aReplicateObj$log2FreqsMatrix, abundanceThresholdsBySampleDf)
# # note: 1 here means sum is applied over rows--i.e., to constructs.
# # Any construct that isn't above sample-specific abundance threshold in at least minNumSamplesAboveAbundance timepoints has
# # 1/TRUE in resulting per-construct vector
# constructFailsMinSamplesAboveAbundanceThresh = apply(passesAbundanceThreshConstructBySampleMatrix, 1,
#                                                      sum) < minNumSamplesAboveAbundance
# # Now amend passing matrix: a construct that passes threshold in fewer than min num required samples should be marked
# # as failing in all samples
# passesAbundanceThreshConstructBySampleMatrix[constructFailsMinSamplesAboveAbundanceThresh,] = FALSE
#
# return (addAbundanceInfoToReplicateObj)
# aReplicateObj$
# # use new matrix to determine which constructs fail abundance threshold across all samples in replicate
#
# aReplicateObj$passesAbundanceThreshMatrix
# # set new matrix as aReplicateData$passesAbundanceThreshMatrix
# # add this new column to aReplicateData$constructAnnotationDf
# }
#
# addAbundanceInfoToReplicateObj < - function(aReplicateObj, passesAbundanceThreshConstructBySampleMatrix,
#                                             constructFailsMinSamplesAboveAbundanceThresh)
# {
# aReplicateObj$passesAbundanceThreshsMatrix = passesAbundanceThreshConstructBySampleMatrix
# aReplicateObj$constructInfoDf = merge(
#     aReplicateObj$constructInfoDf, constructFailsMinSamplesAboveAbundanceThresh, by = "row.names")
# return (aReplicateObj)
# }