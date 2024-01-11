import os
import pandas as pd


def col_from_df(df, column_name) -> list:
    """
    Returns a given dataframe's column as a list object.
    """

    if column_name in df:
        return df[column_name].tolist()
    return []


def csv_to_df(cmd_pointer, filename):
    """
    Returns a dataframe from a csv file.
    """

    if not os.path.isfile(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + filename):
        raise Exception("File does not exist")  # pylint: disable=broad-exception-raised
    else:
        df = pd.read_csv(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + filename)
    return df
