import pandas as pd
from pathlib import Path

class DataLoader():
    """
    Handles loading data from files
    """

    def __init__(self):
        self.data_path = "."

    def set_data_path(self, data_path: str):
        """
        Sets path to get data from
        """
        self.data_path = data_path

    def get_data_path(self):
        """
        Gets data path
        """
        return self.data_path

    def load_files_by_expression(self, expression: str):
        """
        Loads files into dataframes by expression/glob pattern
        """
        files = [file for file in Path(self.get_data_path()).rglob(f"*{expression}*")]
        dfs = [self.load_file(str(file)) for file in files]
        return dfs

    def load_file_by_expression(self, expression: str):
        """
        Loads a file into a dataframe by expression/glob pattern"
        """
        file = [file for file in Path(self.get_data_path()).rglob(f"*{expression}*")][0]
        df = self.load_file(str(file))
        return df

    def load_file(self, filename: str):
        """
        Loads a CSV or XLSX file into a dataframe
        """

        glob: bool = '/' not in filename

        path = Path(filename) if glob == False else [p for p in Path(self.get_data_path()).rglob(f"*{filename}*")][0]
        filepath = path.absolute()

        df: pd.DataFrame = None

        if path.suffix == '.xlsx':
            df = pd.read_excel(pd.ExcelFile(filepath))
        elif path.suffix == '.csv':
            df = pd.read_csv(filepath)
        else:
            raise TypeError("Only XLSX or CSV files supported")

        if "set" in path.name:
            df.name = path.name.split("_")[0]
        else:
            df.name = path.stem + "_" + str(path).split("/")[-2].split("_")[-2] + "_" + str(path).split("/")[-2].split("_")[-1][0]

        return df