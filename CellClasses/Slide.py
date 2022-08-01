import time
from enum import Enum
from typing import List
import os

from CellClasses.Tissue import Tissue


class Slide:

    def __init__(self, data_path: str, genes: List[Enum]):
        self.tissues = None
        self.thresholds = None
        self.data_path = data_path
        self.genes = genes
        self.initialize_tissues()

    def initialize_tissues(self) -> None:
        """
        Initializes the tissues in the slide
        :return: None
        """
        tissues = dict()

        # For each folder in the slide's folder, create a Tissue object
        for folder in [f for f in os.listdir(self.data_path) if not f.startswith('.')]:
            tissues[folder] = Tissue(slide=self, name=folder, data_path=self.data_path + folder)
        self.tissues = tissues

    def get_genes(self):
        """
        Gets the genes for the slide
        :return: list of GeneType
        """
        return self.genes
