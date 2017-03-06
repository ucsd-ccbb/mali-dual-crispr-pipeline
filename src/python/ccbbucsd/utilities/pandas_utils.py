"""This module exposes utility functions and classes for working with pandas objects."""

# third-party libraries
import pandas

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def add_series_to_dataframe(dataframe, series, header):
    """Insert the input series into the input dataframe with the specified column header.

    Args:
        dataframe (pandas.DataFrame): The dataframe to which to add a column; insert is done in-place.
        series (array-like, dict, or scalar value): The column values to add to the dataframe.
        header (str): The name to be used for the new column.
    """
    dataframe.loc[:, header] = pandas.Series(series, index=dataframe.index)


def merge_files_by_shared_header(file_fps, merge_col_header):
    combined_df = None

    for curr_file_fp in file_fps:
        curr_df = pandas.read_table(curr_file_fp)
        if combined_df is None:
            combined_df = curr_df
        else:
            # check to make sure header lists are identical in both files
            if curr_df[[merge_col_header]] != combined_df[[merge_col_header]]:
                raise ValueError("{0} column of file {1} does not match expected list.".format(
                    merge_col_header, curr_file_fp))

            combined_df.merge(curr_df, on=merge_col_header)

        return combined_df
