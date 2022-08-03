from enum import Enum


class RegionType(Enum):
    EverythingCells = "EverythingCells"
    Cytoplasm = "EverythingCytoplasm"
    Nuclei = "EverythingNuclei"
    Extracellular = "EverythingExtracellular"
