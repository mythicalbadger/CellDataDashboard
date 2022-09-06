import pandas as pd
import numpy as np
import itertools
from typing import List

from AuxilaryClasses.DataLoader import DataLoader
import CellClasses.SlideDeck
import CellClasses.RegionType
from CellClasses.RegionType import RegionType
from CellClasses.GeneType import GeneType
#from AuxilaryClasses.DataPlot import DataPlot

import itertools

class Tissue:

    def __init__(self, slide, name: str, data_path: str):
        self.slide = slide
        self.name = name
        self.data_path = data_path
        self.datasets = dict()
        self.data_loader = DataLoader()
        self.data_loader.set_data_path(data_path)
        self.load_data()
        self.db = None

        self.everything_cells = self.sort_datasets_by_field(self.get_datasets_by_pattern("EverythingCells_set") + self.get_datasets_by_pattern("EverythingCells_Field"))
        self.cytoplasm = self.sort_datasets_by_field(self.get_datasets_by_pattern("EverythingCytoplasm_set") + self.get_datasets_by_pattern("EverythingCytoplasm_Field"))
        self.nuclei = self.sort_datasets_by_field(self.get_datasets_by_pattern("EverythingNuclei_set") + self.get_datasets_by_pattern("EverythingNuclei_Field"))
        self.extracellular = self.sort_datasets_by_field(self.get_datasets_by_pattern("EverythingExtracellular_set") + self.get_datasets_by_pattern("EverythingExtracellular_Field"))
        self.region_map = {
            RegionType.EverythingCells: self.everything_cells,
            RegionType.Cytoplasm: self.cytoplasm,
            RegionType.Nuclei: self.nuclei,
            RegionType.Extracellular: self.extracellular
        }

        self.initialize_db()

    def db_contains(self, field: str, region: RegionType):
        return self.db[region][field] is not None

    def db_put(self, field: str, data: pd.DataFrame, region: RegionType):
        self.db[region][field] = data

    def db_get(self, field: str, region: RegionType):
        return self.db[region][field]

    def initialize_db(self):
        fields = [
            "CellCounts",
            "MedianMarkers",
            "ExpressedMarkers",
            "ACDScores",
            "ZeroScores",
            "PositiveExpressors",
            "PercentBins",
            "WeightedPercentBins",
            "HScores",
            "HScoreMedians",
            "BorderHScoreMedians",
            "MiddleHScoreMedians",
            "MedianNbrs",
            "MedianPercentTouching",
            "PairGeneExpression",
            "TrioGeneExpression",
            "QuadGeneExpression"
        ]

        db = dict()
        for region in self.region_map:
            db[region] = {f: None for f in fields}
        self.db = db

    def load_data(self) -> None:
        """
        Loads required datasets -- EverythingExpression and EverythingCells for each field
        :return: None
        """
        self.add_datasets(self.data_loader.load_files_by_expression("EverythingExpression"))
        self.add_datasets(self.data_loader.load_files_by_expression("EverythingCells"))
        self.add_datasets(self.data_loader.load_files_by_expression("EverythingNuclei"))
        self.add_datasets(self.data_loader.load_files_by_expression("EverythingCytoplasm"))
        self.add_datasets(self.data_loader.load_files_by_expression("EverythingExtracellular"))

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

    def filter_datasets_by_name_pattern(self, datasets: List[pd.DataFrame], pattern: str):
        return list(filter(lambda d: pattern in d.name, datasets))

    def filter_datasets_by_num_token(self, datasets: List[pd.DataFrame], token: str, num: int):
        return list(filter(lambda d: d.name.count(token) == num, datasets))

    def get_pair_gene_datasets(self):
        return self.filter_datasets_by_num_token(list(self.datasets.values()), "and", 1)

    def get_trio_gene_datasets(self):
        return self.filter_datasets_by_num_token(list(self.datasets.values()), "and", 2)

    def get_quad_gene_datasets(self):
        return self.filter_datasets_by_num_token(list(self.datasets.values()), "and", 3)

    def list_datasets(self) -> None:
        """
        Lists out all datasets for this tissue
        :return: None
        """
        for i, k in enumerate(self.datasets):
            print(f"{i + 1}) {k}")

    def sort_datasets_by_field(self, datasets):
        """
        Sorts datasets by their field number
        :param datasets: a list of datasets
        :return: a list of datasets sorted by field number
        """
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

    def region_to_data(self, region: RegionType):
        return self.region_map[region]

    def calculate_cell_counts(self, region: RegionType):
        """
        Calculates total cell counts by fetching EverythingCells
        :return: pandas DataFrame with counts
        """
        if self.db_contains("CellCounts", region):
            return self.db_get("CellCounts", region)

        datasets = self.region_to_data(region)

        # Name format: Exp_EverythingCells_Field <1-14>_<B/M>
        fields = [d.name.split(f"{region.value}_")[-1] for d in datasets]  # field names
        lengths = [len(d.index) for d in datasets]  # cell counts

        ret = pd.DataFrame([lengths], index=["Counts"], columns=fields)
        self.db_put("CellCounts", ret, region)
        return ret

    def calculate_expressed_markers(self, region: RegionType):
        """
        Calculates total expressed markers by fetching EverythingCells
        :return: DataFrame of total expressed markers
        """
        if self.db_contains("ExpressedMarkers", region):
            return self.db_get("ExpressedMarkers", region)

        # Get dataset and column key from file with all cell data (nuclei + cytoplasm)
        datasets = self.region_to_data(region)

        # Generate multiindex with gene name on the outside and ACD score field on the inside
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["Expressions"])
        scores = {d.name.split(f"{region.value}_")[-1]: [] for d in datasets}

        # Calculate ACD scores for each gene for each dataset
        for gene in self.slide.get_genes():
            for dataset in datasets:
                column_key = "Children_Expression_%s_Count" % gene.value
                scores[dataset.name.split(f"{region.value}_")[-1]].append(dataset.loc[:, column_key].sum())

        ret = pd.DataFrame(scores, index=idx)
        self.db_put("ExpressedMarkers", ret, region)
        return ret

    def calculate_median_markers(self, region: RegionType):
        """
        Calculates the median expressed markers for each cell (by gene)
        :return: DataFrame of median expressed markers
        """
        if self.db_contains("MedianMarkers", region):
            return self.db_get("MedianMarkers", region)

        # Get dataset and column key from file with all cell data (nuclei + cytoplasm)
        datasets = self.region_to_data(region)

        # Generate multiindex with gene name on the outside and ACD score field on the inside
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["Expressions"])
        scores = {d.name.split(f"{region.value}_")[-1]: [] for d in datasets}

        # Calculate ACD scores for each gene for each dataset
        for gene in self.slide.get_genes():
            for dataset in datasets:
                column_key = "Children_Expression_%s_Count" % gene.value
                scores[dataset.name.split(f"{region.value}_")[-1]].append(dataset.loc[:, column_key].median())

        ret = pd.DataFrame(scores, index=idx)
        self.db_put("MedianMarkers", ret, region)
        return ret

    def score_acd_ranges(self, region: RegionType):
        """
        Scores ACD ranges for each field
        - ACD Score 0 : number of cells with no expressions
        - ACD Score 1 : number of cells with expressions equal to the threshold for the corresponding gene
        - ACD Score 2 : number of cells with expressions greater than ACD Score 1 threshold and less or equal to nine
        - ACD Score 3 : number of cells with expressions greater than or equal to 10 and less than or equal to 15
        - ACD Score 4 : number of cells with expressions greater than 15
        :return: pandas DataFrame with ACD scores
        """
        if self.db_contains("ACDScores", region):
            return self.db_get("ACDScores", region)

        # Get dataset and column key from file with all cell data (nuclei + cytoplasm)
        datasets = self.region_to_data(region)

        # Generate multiindex with gene name on the outside and ACD score field on the inside
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["ACD Score 1", "ACD Score 2", "ACD Score 3", "ACD Score 4"])
        scores = {d.name.split(f"{region.value}_")[-1]: [] for d in datasets}

        # Calculate ACD scores for each gene for each dataset
        for gene in self.slide.get_genes():
            for dataset in datasets:
                column_key = "Children_Expression_%s_Count" % gene.value

                # Gets expression counts
                counts = dataset[column_key].value_counts().sort_index()
                indices = counts.index.tolist()

                # Merges counts and indices into a pandas dataframe
                count_df = pd.DataFrame(zip(indices, counts), columns=["Expressions", "Score"])

                # Calculate range scores
                range_scores = [
                    sum(count_df.loc[(count_df['Expressions'] >= 1) & (count_df['Expressions'] <= 3)]['Score']),
                    sum(count_df.loc[(count_df['Expressions'] >= 4) & (count_df['Expressions'] <= 9)]['Score']),
                    sum(count_df.loc[(count_df['Expressions'] >= 10) & (count_df['Expressions'] <= 15)]['Score']),
                    sum(count_df.loc[count_df['Expressions'] > 15]['Score']),
                ]

                scores[dataset.name.split(f"{region.value}_")[-1]].extend(range_scores)

        ret = pd.DataFrame(scores, index=idx)
        self.db_put("ACDScores", ret, region)
        return ret

    def calculate_zero_scores(self, region: RegionType):
        """
        Calculates zero scores for each field - the total number of expressions subtracted from the total number of cells
        :return: pandas DataFrame with zero scores
        """
        if self.db_contains("ZeroScores", region):
            return self.db_get("ZeroScores", region)

        # Get dataset and column key from file with all cell data (nuclei + cytoplasm)
        datasets = self.region_to_data(region)

        # Generate multiindex with gene name on the outside and ACD score field on the inside
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["ACD Score 0"])

        scores = {d.name.split(f"{region.value}_")[-1]: [] for d in datasets}

        # Calculate ACD scores for each gene for each dataset
        for i, dataset in enumerate(datasets):
            for gene in self.slide.get_genes():
                column_key = "Children_Expression_%s_Count" % gene.value

                # Gets expression counts
                zsc = dataset[column_key].where(lambda x: x == 0).dropna().size
                scores[dataset.name.split(f"{region.value}_")[-1]].append(zsc)

        ret = pd.DataFrame(scores, index=idx)
        self.db_put("ZeroScores", ret, region)
        return ret

    def calculate_positive_expressors(self, region: RegionType):
        """
        Calculates the number of cells that have expressions >= bin_one_threshold
        :return:
        """
        if self.db_contains("PositiveExpressors", region):
            return self.db_get("PositiveExpressors", region)

        # Get dataset and column key from file with all cell data (nuclei + cytoplasm)
        datasets = self.region_to_data(region)

        # Generate multiindex with gene name on the outside and ACD score field on the inside
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["Positive Expressors"])
        scores = {d.name.split(f"{region.value}_")[-1]: [] for d in datasets}

        # Calculate ACD scores for each gene for each dataset
        for i, dataset in enumerate(datasets):
            for gene in self.slide.get_genes():
                column_key = "Children_Expression_%s_Count" % gene.value

                # Gets expression counts
                bin_one_threshold = CellClasses.SlideDeck.SlideDeck.get_threshold(gene)
                zsc = dataset[column_key].where(lambda x: x >= bin_one_threshold).dropna().size
                scores[dataset.name.split(f"{region.value}_")[-1]].append(zsc)

        ret = pd.DataFrame(scores, index=idx)
        self.db_put("PositiveExpressors", ret, region)
        return ret

    def score_percent_bins(self, region: RegionType):
        """
        Calculates percent bins for each field
        - Percent bin 0 : percentage of number of cells with 0 expressions in the dataset
        - Percent bin 1 : percentage of number of cells with ACD score 1 in the dataset
        - Percent bin 2 : percentage of number of cells with ACD score 2 in the dataset
        - Percent bin 3 : percentage of number of cells with ACD score 3 in the dataset
        - Percent bin 4 : percentage of number of cells with ACD score 4 in the dataset
        :return: pandas DataFrame with percent bin scores
        """
        if self.db_contains("PercentBins", region):
            return self.db_get("PercentBins", region)

        # Get pertinent datasets and setup return dataframe
        datasets = self.region_to_data(region)
        acd_scores = self.db_get("ACDScores", region) if self.db_contains("ACDScores", region) else self.score_acd_ranges(region)
        zero_scores = self.db_get("ZeroScores", region) if self.db_contains("ZeroScores", region) else self.calculate_zero_scores(region)
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

        ret = pd.DataFrame(percents, index=idx)
        self.db_put("PercentBins", ret, region)
        return ret

    def score_weighted_percent_bins(self, region: RegionType):
        """
        Calculates weighted percent bins for each field. Basically just (percent bin score X bin number)
        :return: pandas DataFrame with weighted percent bin scores
        """
        if self.db_contains("WeightedPercentBins", region):
            return self.db_get("WeightedPercentBins", region)

        # Get percent bins and set up multipliers / index
        percent_bins = self.db_get("PercentBins", region) if self.db_contains("PercentBins", region) else self.score_percent_bins(region)
        mult = [i for i in range(5)] * len(self.slide.get_genes())
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], [f"Weighted bin {i}" for i in range(5)])
        weighted = {col: [] for col in percent_bins.columns}

        # Multiply each column by mult to get weighted
        for i, col in enumerate(percent_bins.columns):
            percent_bins[col] *= mult
            weighted[col] = percent_bins.loc[:, col].tolist()

        ret = pd.DataFrame(weighted, index=idx)
        self.db_put("WeightedPercentBins", ret, region)
        return ret

    def calculate_hscores(self, region: RegionType):
        """
        Calculates H-Scores for each gene - the sum of weighted percent bins
        :return: pandas DataFrame with H-Scores
        """
        if self.db_contains("HScores", region):
            return self.db_get("HScores", region)

        # Get weighted percent bins and set up index
        weighted = self.db_get("WeightedPercentBins", region) if self.db_contains("WeightedPercentBins", region) else self.score_weighted_percent_bins(region)
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["H-Score"])
        hscores = {col: [] for col in weighted.columns}

        for i, col in enumerate(weighted.columns):
            field = weighted[col].tolist()
            hscores[col] = [sum(field[i:i + 5]) for i in range(0, len(field), 5)]

        ret = pd.DataFrame(hscores, index=idx)
        self.db_put("HScores", ret, region)
        return ret

    def calculate_hscore_medians(self, region: RegionType):
        """
        Calculates the total median of H-Scores
        :return: pandas DataFrame containing median
        """
        if self.db_contains("HScoreMedians", region):
            return self.db_get("HScoreMedians", region)

        hscores = self.db_get("HScores", region) if self.db_contains("HScores", region) else self.calculate_hscores(region)
        ret = hscores.median(axis=1)
        self.db_put("HScoreMedians", ret, region)
        return ret

    def calculate_border_hscore_medians(self, region: RegionType):
        """
        Calculates the border medians of H-Scores (fields with B)
        :return: pandas DataFrame containing medians
        """
        if self.db_contains("BorderHScoreMedians", region):
            return self.db_get("BorderHScoreMedians", region)
        hscores = self.db_get("HScores", region) if self.db_contains("HScores", region) else self.calculate_hscores(region)
        ret = hscores.loc[:, hscores.columns.str.contains("_B")].median(axis=1)
        self.db_put("BorderHScoreMedians", ret, region)
        return ret

    def calculate_middle_hscore_medians(self, region: RegionType):
        """
        Calculates the middle median of H-Scores (fields with M)
        :return: pandas DataFrame containing medians
        """
        if self.db_contains("MiddleHScoreMedians", region):
            return self.db_get("MiddleHScoreMedians", region)
        hscores = self.db_get("HScores", region) if self.db_contains("HScores", region) else self.calculate_hscores(region)
        ret = hscores.loc[:, hscores.columns.str.contains("_M")].median(axis=1)
        self.db_put("MiddleHScoreMedians", ret, region)
        return ret

    def hscores_to_csv(self, region):
        total = self.calculate_hscore_medians(region)
        total_idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["Total H-Scores"])
        total.index = total_idx

        border = self.calculate_border_hscore_medians(region)
        border_idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], ["Border H-Scores"])
        border.index = border_idx

        middle = self.calculate_middle_hscore_medians(region)
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

    def median_nbrs(self, region: RegionType):
        """
        Calculates the median number of neighbors of each gene type per gene per cell(???)
        """
        if self.db_contains("MedianNbrs", region):
            return self.db_get("MedianNbrs", region)

        gene_to_num = {g.value : i for i, g in enumerate(self.slide.get_genes())}
        if GeneType.NCBP2AS2 in self.slide.get_genes():
            gene_to_num["NCBP2"] = gene_to_num.pop(GeneType.NCBP2AS2.value)
        # I hate myself
        datasets = sorted(self.sort_datasets_by_field(list(filter(lambda d: f"{region.value}_Field" not in d.name and "and" not in d.name, self.get_datasets_by_pattern(f"Exp_{region.value}")))), key=lambda d: gene_to_num[d.name.split("_")[2]])
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], [g.value for g in self.slide.get_genes()])
        medians = {"Field " + d.name.split(" ")[-1]: [] for d in datasets}  # field names
        for dataset in datasets:
            for gene in self.slide.get_genes():
                column_key = f"Neighbors_NumberOfNeighbors_Expression_{gene.value}_20"
                median = 0.0
                if column_key in dataset.columns:
                    median = dataset.loc[:, column_key].median() if not dataset.loc[:, column_key].isna().all() else 0.0
                medians["Field " + dataset.name.split(" ")[-1]].append(median)
        ret = pd.DataFrame(medians, index=idx)
        self.db_put("MedianNbrs", ret, region)
        return ret

    def median_percent_touching(self, region: RegionType):
        """
        Calculates the median number of neighbors touching for each gene type per gene per cell(???)
        """
        if self.db_contains("MedianPercentTouching", region):
            return self.db_get("MedianPercentTouching", region)

        gene_to_num = {g.value : i for i, g in enumerate(self.slide.get_genes())}
        if GeneType.NCBP2AS2 in self.slide.get_genes():
            gene_to_num["NCBP2"] = gene_to_num.pop(GeneType.NCBP2AS2.value)
        # I hate myself
        datasets = sorted(self.sort_datasets_by_field(list(filter(lambda d: f"{region.value}_Field" not in d.name and "and" not in d.name, self.get_datasets_by_pattern(f"Exp_{region.value}")))), key=lambda d: gene_to_num[d.name.split("_")[2]])
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], [g.value for g in self.slide.get_genes()])
        percents = {"Field " + d.name.split(" ")[-1]: [] for d in datasets}  # field names
        for dataset in datasets:
            for gene in self.slide.get_genes():
                column_key = f"Neighbors_PercentTouching_Expression_{gene.value}_20"
                percent = 0.0
                if column_key in dataset.columns:
                    percent = dataset.loc[:, column_key].median() if not dataset.loc[:, column_key].isna().all() else 0.0
                percents["Field " + dataset.name.split(" ")[-1]].append(percent)
        ret = pd.DataFrame(percents, index=idx)
        self.db_put("PercentTouching", ret, region)
        return ret


    def calculate_pair_gene_expression(self, region: RegionType):
        """
        Calculates gene pair counts by fetching EverythingCells
        :return: pandas DataFrame with counts
        """
        if self.db_contains("PairGeneExpression", region):
            return self.db_get("PairGeneExpression", region)

        datasets = self.region_to_data(region)
        expressions = {d.name.split(f"{region.value}_")[-1]: [] for d in datasets}

        pair_cols = list(filter(lambda col_name: col_name.count("and") == 1, datasets[0].columns.tolist()))
        if len(pair_cols) == 0:
            return pd.DataFrame([])

        idx = self.gen_multi_idx([col[15:-6] for col in pair_cols], ["Expressions"])

        for dataset in datasets:
            for col in pair_cols:
                s = dataset.loc[:, col].where(lambda x: x > 0).dropna().count()
                expressions[dataset.name.split(f"{region.value}_")[-1]].append(s)

        ret = pd.DataFrame(expressions, index=idx)
        self.db_put("PairGeneExpression", ret, region)
        return ret


    def calculate_trio_gene_expression(self, region: RegionType):
        """
        Calculates gene trio counts by fetching EverythingCells
        :return: pandas DataFrame with counts
        """
        if self.db_contains("TrioGeneExpression", region):
            return self.db_get("TrioGeneExpression", region)

        datasets = self.region_to_data(region)
        expressions = {d.name.split(f"{region.value}_")[-1]: [] for d in datasets}

        pair_cols = list(filter(lambda col_name: col_name.count("and") == 2, datasets[0].columns.tolist()))
        if len(pair_cols) == 0:
            return pd.DataFrame([])

        idx = self.gen_multi_idx([col[15:-6] for col in pair_cols], ["Expressions"])

        for dataset in datasets:
            for col in pair_cols:
                s = dataset.loc[:, col].where(lambda x: x > 0).dropna().count()
                expressions[dataset.name.split(f"{region.value}_")[-1]].append(s)

        ret = pd.DataFrame(expressions, index=idx)
        self.db_put("TrioGeneExpression", ret, region)
        return ret

    def calculate_quad_gene_expression(self, region: RegionType):
        """
        Calculates gene quadruple counts by fetching EverythingCells
        :return: pandas DataFrame with counts
        """
        if self.db_contains("QuadGeneExpression", region):
            return self.db_get("QuadGeneExpression", region)

        datasets = self.region_to_data(region)
        expressions = {d.name.split(f"{region.value}_")[-1]: [] for d in datasets}

        pair_cols = list(filter(lambda col_name: col_name.count("and") == 3, datasets[0].columns.tolist()))
        if len(pair_cols) == 0:
            return pd.DataFrame([])

        idx = self.gen_multi_idx([col.split("_")[-2] for col in pair_cols], ["Expressions"])

        for dataset in datasets:
            for col in pair_cols:
                s = dataset.loc[:, col].where(lambda x: x > 0).dropna().count()
                expressions[dataset.name.split(f"{region.value}_")[-1]].append(s)

        ret = pd.DataFrame(expressions, index=idx)
        self.db_put("QuadGeneExpression", ret, region)
        return ret

    def euclidean_distance(self, p1, p2):
        xs1, ys1 = p1
        xs2, ys2 = p2
        return np.sqrt((xs1-xs2)**2 + (ys1-ys2)**2)

    def median_distance(self, dataA, dataB):
        distances = [self.euclidean_distance(A, B) for A in dataA for B in dataB]
        return np.median(np.array(distances))

    """
    def calculate_median_distances(self, region: RegionType):
        datasets = self.region_to_data(region)
        ret = { d.name.split(f"{region.value}_")[-1]: [] for d in datasets}
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], [g.value for g in self.slide.get_genes()])

        for df in datasets:
            plots = []
            for gene in self.slide.get_genes():
                gene_dataset = df.name.split("_")
                gene_dataset.insert(-2, gene.value)
                gene_dataset = '_'.join(gene_dataset)
                plot = DataPlot(
                    gene.value,
                    np.array(self.get_dataset(gene_dataset).loc[:, f"Location_Center_X"].dropna()), 
                    np.array(self.get_dataset(gene_dataset).loc[:, f"Location_Center_Y"].dropna()), 
                )
                plots.append(plot)

            for i in range(len(plots)):
                for j in range(len(plots)):
                    g1 = plots[i]
                    g2 = plots[j]
                    g1_data = list(zip(g1.xs, g1.ys))                
                    g2_data = list(zip(g2.xs, g2.ys))                
                    ret[df.name.split(f"{region.value}_")[-1]].append(self.median_distance(g1_data, g2_data))
            
        ret = pd.DataFrame(ret, index=idx)
        return ret

    def calculate_median_pair_distances(self, region: RegionType):
        datasets = self.region_to_data(region)
        ret = { d.name.split(f"{region.value}_")[-1]: [] for d in datasets}
        pair_datasets = self.get_pair_gene_datasets()
        print([df.name for df in pair_datasets])
        idx = self.gen_multi_idx([g.value for g in self.slide.get_genes()], [g.value for g in self.slide.get_genes()])

        for df in datasets:
            plots = []
            for gene in self.slide.get_genes():
                gene_dataset = df.name.split("_")
                gene_dataset.insert(-2, gene.value)
                gene_dataset = '_'.join(gene_dataset)
                plot = DataPlot(
                    gene.value,
                    np.array(self.get_dataset(gene_dataset).loc[:, f"Location_Center_X"].dropna()), 
                    np.array(self.get_dataset(gene_dataset).loc[:, f"Location_Center_Y"].dropna()), 
                )
                plots.append(plot)

            for i in range(len(plots)):
                for j in range(len(plots)):
                    g1 = plots[i]
                    g2 = plots[j]
                    g1_data = list(zip(g1.xs, g1.ys))                
                    g2_data = list(zip(g2.xs, g2.ys))                
                    ret[df.name.split(f"{region.value}_")[-1]].append(self.median_distance(g1_data, g2_data))
            
        ret = pd.DataFrame(ret, index=idx)
        return ret
        """
