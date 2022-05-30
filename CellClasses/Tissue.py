import pandas as pd
from typing import List
import numpy as np

from AuxilaryClasses.DataLoader import DataLoader
import CellClasses.SlideDeck


class Tissue:

    def __init__(self, slide, name: str, data_path: str):
        self.slide = slide
        self.name = name
        self.data_path = data_path
        self.datasets = dict()
        self.data_loader = DataLoader()
        self.data_loader.set_data_path(data_path)
        self.load_data()

    def load_data(self) -> None:
        """
        Loads required datasets -- EverythingExpression and EverythingCells for each field
        :return: None
        """
        self.add_datasets(self.data_loader.load_files_by_expression("EverythingExpression"))
        self.add_datasets(self.data_loader.load_files_by_expression("EverythingCells.csv"))

    def add_dataset(self, dataset: pd.DataFrame) -> None:
        """
        Adds an individual dataset to this tissue's list of datasets
        :param dataset: the pandas DataFrame to add
        :return: None
        """
        self.datasets[dataset.name] = dataset

    def add_datasets(self, datasets: List[pd.DataFrame]) -> None:
        """
        Adds a list of datasets to this tissue's lis of datasets
        :param datasets: the list of pandas DataFrames to add
        :return: None
        """
        for dataset in datasets:
            self.add_dataset(dataset)

    def get_dataset(self, dataset_name: str) -> pd.DataFrame:
        """
        Fetches a dataset by name
        :param dataset_name: the name of the dataset to fetch
        :return: the matching pandas DataFrame
        """
        return self.datasets[dataset_name]

    def get_datasets_by_pattern(self, pattern: str) -> List[pd.DataFrame]:
        """
        Returns a list of datasets whose names match a specific pattern
        :param pattern: the pattern to match
        :return: a list of pandas DataFrames
        """
        return [self.get_dataset(d) for d in self.datasets if pattern in d]

    def list_datasets(self) -> None:
        """
        Lists out all datasets for this tissue
        :return: None
        """
        for i, k in enumerate(self.datasets):
            print(f"{i + 1}) {k}")

    def sort_datasets_by_field(self, datasets):
        return sorted(datasets, key=lambda d: int(d.name.split(" ")[-1].split("_")[0]))

    def gen_multi_idx(self, left, right):
        """
        Generates a pandas multi-index for a dataframe
        :param left: the leftmost, outer index that covers multiple rows
        :param right: the inner index that is in every row
        :return: pandas MultiIndex object
        """
        indices = [(l, r) for l in left for r in right]
        return pd.MultiIndex.from_tuples(indices)

    def calculate_cell_counts(self):
        """
        Calculates total cell counts by fetching EverythingCells
        :return: pandas DataFrame with counts
        """

        # Get EverythingCells datasets for each field
        datasets = self.sort_datasets_by_field(self.get_datasets_by_pattern("EverythingCells"))

        # Name format: Exp_EverythingCells_Field <1-14>_<B/M>
        # Just want the Field <1-14>
        fields = [d.name.split("Cells_")[-1] for d in datasets]

        # Number of rows are the counts
        lengths = [len(d.index) for d in datasets]

        # Zip them together and sort by field number
        counts = np.array(list(zip(fields, lengths)))

        # Turn them into dictionary for pandas
        df = {counts[:, 0][i]: counts[:, 1][i] for i in range(len(counts))}

        return pd.DataFrame(df, index=["Counts"])

    def calculate_expressed_markers(self):
        """
        Calculates total expressed markers by fetching EverythingCells
        """
        # Get dataset and column key from file with all cell data (nuclei + cytoplasm)
        datasets = self.sort_datasets_by_field(self.get_datasets_by_pattern("EverythingCells"))

        # Generate multiindex with gene name on the outside and ACD score field on the inside
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()],
                                 ["Expressions"])

        scores = {datasets[i].name.split("Cells_")[-1]: [] for i in range(len(datasets))}

        # Calculate ACD scores for each gene for each dataset
        for gene in self.slide.get_genes():
            for i, dataset in enumerate(datasets):
                column_key = "Children_Expression_%s_Count" % gene.value
                scores[dataset.name.split("Cells_")[-1]].append(dataset.loc[:, column_key].sum())

        ret = pd.DataFrame(scores, index=idx)
        return ret

    def score_acd_ranges(self):
        """
        Scores ACD ranges for each field
        - ACD Score 0 : number of cells with no expressions
        - ACD Score 1 : number of cells with expressions equal to the threshold for the corresponding gene
        - ACD Score 2 : number of cells with expressions greater than ACD Score 1 threshold and less or equal to nine
        - ACD Score 3 : number of cells with expressions greater than or equal to 10 and less than or equal to 15
        - ACD Score 4 : number of cells with expressions greater than 15
        :return: pandas DataFrame with ACD scores
        """
        # Get dataset and column key from file with all cell data (nuclei + cytoplasm)
        datasets = self.sort_datasets_by_field(self.get_datasets_by_pattern("EverythingCells"))

        # Generate multiindex with gene name on the outside and ACD score field on the inside
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()],
                                 ["ACD Score 1", "ACD Score 2", "ACD Score 3", "ACD Score 4"])

        scores = {datasets[i].name.split("Cells_")[-1]: [] for i in range(len(datasets))}

        # Calculate ACD scores for each gene for each dataset
        for gene in self.slide.get_genes():
            for i, dataset in enumerate(datasets):
                column_key = "Children_Expression_%s_Count" % gene.value

                # Gets expression counts
                counts = dataset[column_key].value_counts().sort_index()
                indices = counts.index.tolist()

                # Merges counts and indices into a pandas dataframe
                count_df = pd.DataFrame(zip(indices, counts), columns=["Expressions", "Score"])

                # Calculate range scores
                range_scores = [
                    sum(count_df.loc[(count_df['Expressions'] >= 1) & (count_df['Expressions'] <= 3)]['Score']),
                    sum(count_df.loc[
                            (count_df['Expressions'] >= 4) & (count_df['Expressions'] <= 9)][
                            'Score']),
                    sum(count_df.loc[(count_df['Expressions'] >= 10) & (count_df['Expressions'] <= 15)]['Score']),
                    sum(count_df.loc[count_df['Expressions'] > 15]['Score']),
                ]

                scores[dataset.name.split("Cells_")[-1]].extend(range_scores)

        ret = pd.DataFrame(scores, index=idx)
        return ret

    def calculate_zero_scores(self):
        """
        Calculates zero scores for each field - the total number of expressions subtracted from the total number of cells
        :return: pandas DataFrame with zero scores
        """
        # Get dataset and column key from file with all cell data (nuclei + cytoplasm)
        datasets = self.sort_datasets_by_field(self.get_datasets_by_pattern("EverythingCells"))

        # Generate multiindex with gene name on the outside and ACD score field on the inside
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()],
                                 ["ACD Score 0"])

        scores = {datasets[i].name.split("Cells_")[-1]: [] for i in range(len(datasets))}

        # Calculate ACD scores for each gene for each dataset
        for i, dataset in enumerate(datasets):
            for gene in self.slide.get_genes():
                column_key = "Children_Expression_%s_Count" % gene.value

                # Gets expression counts
                zsc = dataset[column_key].where(lambda x: x == 0).dropna().size
                scores[dataset.name.split("Cells_")[-1]].append(zsc)

        ret = pd.DataFrame(scores, index=idx)
        return ret

    def calculate_positive_expressors(self):
        # Get dataset and column key from file with all cell data (nuclei + cytoplasm)
        datasets = self.sort_datasets_by_field(self.get_datasets_by_pattern("EverythingCells"))

        # Generate multiindex with gene name on the outside and ACD score field on the inside
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()],
                                 ["Positive Expressors"])

        scores = {datasets[i].name.split("Cells_")[-1]: [] for i in range(len(datasets))}

        # Calculate ACD scores for each gene for each dataset
        for i, dataset in enumerate(datasets):
            for gene in self.slide.get_genes():
                column_key = "Children_Expression_%s_Count" % gene.value

                # Gets expression counts
                bin_one_threshold = CellClasses.SlideDeck.SlideDeck.get_threshold(gene)
                zsc = dataset[column_key].where(lambda x: x >= bin_one_threshold).dropna().size
                scores[dataset.name.split("Cells_")[-1]].append(zsc)

        ret = pd.DataFrame(scores, index=idx)
        return ret

    def score_percent_bins(self):
        """
        Calculates percent bins for each field
        - Percent bin 0 : percentage of number of cells with 0 expressions in the dataset
        - Percent bin 1 : percentage of number of cells with ACD score 1 in the dataset
        - Percent bin 2 : percentage of number of cells with ACD score 2 in the dataset
        - Percent bin 3 : percentage of number of cells with ACD score 3 in the dataset
        - Percent bin 4 : percentage of number of cells with ACD score 4 in the dataset
        :return: pandas DataFrame with percent bin scores
        """

        # Get pertinent datasets and setup return dataframe
        datasets = self.sort_datasets_by_field(self.get_datasets_by_pattern("EverythingCells"))
        acd_scores = self.score_acd_ranges()
        zero_scores = self.calculate_zero_scores()
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], [f"% bin {i}" for i in range(5)])
        percents = {col: [] for col in acd_scores.columns}

        # Calculate the percent bins for ACD Score 1 - 4
        for i, col in enumerate(acd_scores.columns):
            percent = acd_scores.loc[:, col] / len(datasets[i]) * 100
            percents[col].extend(percent)

        # Calculate the percent bins for ACD Score 0
        for i, percent in enumerate(percents):
            for j in range(0, len(acd_scores), 5):
                zs = zero_scores.loc[:, acd_scores.columns[i]][0] / len(datasets[i]) * 100
                percents[percent].insert(j, zs)

        return pd.DataFrame(percents, index=idx)

    def score_weighted_percent_bins(self):
        """
        Calculates weighted percent bins for each field. Basically just (percent bin score X bin number)
        :return: pandas DataFrame with weighted percent bin scores
        """

        # Get percent bins and set up multipliers / index
        percent_bins = self.score_percent_bins()
        mult = [i for i in range(5)] * len(self.slide.get_genes())
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], [f"Weighted bin {i}" for i in range(5)])
        weighted = {col: [] for col in percent_bins.columns}

        # Multiply each column by mult to get weighted
        for i, col in enumerate(percent_bins.columns):
            percent_bins[col] *= mult
            weighted[col] = percent_bins.loc[:, col].tolist()

        return pd.DataFrame(weighted, index=idx)

    def calculate_hscores(self):
        """
        Calculates H-Scores for each gene - the sum of weighted percent bins
        :return: pandas DataFrame with H-Scores
        """

        # Get weighted percent bins and set up index
        weighted = self.score_weighted_percent_bins()
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["H-Score"])
        hscores = {col: [] for col in weighted.columns}

        for i, col in enumerate(weighted.columns):
            field = weighted[col].tolist()
            hscores[col] = [sum(field[i:i + 5]) for i in range(0, len(field), 5)]

        return pd.DataFrame(hscores, index=idx)

    def calculate_hscore_means(self):
        """
        Calculates the total mean of H-Scores
        :return: pandas DataFrame containing means
        """
        hscores = self.calculate_hscores()
        return hscores.mean(axis=1)

    def calculate_border_hscore_means(self):
        """
        Calculates the border mean of H-Scores (fields with B)
        :return: pandas DataFrame containing means
        """
        hscores = self.calculate_hscores()
        return hscores.loc[:, hscores.columns.str.contains("_B")].mean(axis=1)

    def calculate_middle_hscore_means(self):
        """
        Calculates the middle mean of H-Scores (fields with M)
        :return: pandas DataFrame containing means
        """
        hscores = self.calculate_hscores()
        return hscores.loc[:, hscores.columns.str.contains("_M")].mean(axis=1)

    def hscores_to_csv(self):
        total = self.calculate_hscore_means()
        total_idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["Total H-Scores"])
        total.index = total_idx

        border = self.calculate_border_hscore_means()
        border_idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["Border H-Scores"])
        border.index = border_idx

        middle = self.calculate_middle_hscore_means()
        middle_idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["Middle H-Scores"])
        middle.index = middle_idx

        frames = [total, border, middle]
        valid_frames = []

        for frame in frames:
            if frame.isnull().values.any(): continue
            valid_frames.append(frame)
        return pd.concat(valid_frames)

    def means_to_csv(self):
        return self.calculate_hscores()


