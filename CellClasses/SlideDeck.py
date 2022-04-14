import os
import CellClasses.Slide
from CellClasses.GeneType import *
import pandas as pd

from pathlib import Path


class SlideDeck:

    def __init__(self, data_path: str) -> None:
        self.slides = dict()
        self.data_path = data_path
        SlideDeck.thresholds: pd.DataFrame = pd.read_excel(pd.ExcelFile(
            data_path + "Thresholds.xlsx"), 'opal color')

        genes = [gene for gene in GeneType]
        self.slide_to_gene = {
            "QC Slide": genes[0:4],
            "Slide 1": genes[4:8],
            "Slide 2": genes[8:12],
            "Slide 3": genes[12:16],
            "Slide 4": genes[16:20],
            "Slide 5": genes[20:24],
            "Slide 6": genes[24:28]
        }

        self.initialize_slides()

    def initialize_slides(self) -> None:
        """
        Initializes the slides in the deck
        :return: None
        """
        slides = dict()
        for folder in sorted(os.listdir(self.data_path)):
            path: str = str(Path(self.data_path + folder).absolute()) + "/"
            slide = CellClasses.Slide.Slide(path, self.slide_to_gene[folder])
            slides[folder] = slide

        self.slides = slides

    @staticmethod
    def get_threshold(gene: GeneType) -> int:
        """
        Gets threshold value for a specific gene
        :param gene: the GeneType of the gene
        :return: integer threshold value
        """
        thresholds: pd.DataFrame = SlideDeck.thresholds

        # "Unnamed: 2" column has gene name
        gene_row: pd.DataFrame = thresholds.loc[thresholds['Unnamed: 2'] == gene.value]
        return int(gene_row.iloc[:, -1])

    def get_slide(self, slide_name: str):
        """
        Fetches a Slide from the deck
        :param slide_num: the number of the slide to get
        :return: matching Slide object
        """
        return self.slides[slide_name]
