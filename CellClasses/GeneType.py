from enum import Enum

class GeneType(Enum):
    # QC
    HPRT1 = "HPRT1"
    POLR2A = "POLR2A"
    PPIB = "PPIB"
    UBC = "UBC"

    # Slide 1
    COL1A2 = "COL1A2"
    PDGFRA = "PDGFRA"
    FOXP3 = "FOXP3"
    CXCL12 = "CXCL12"

    # Slide 2
    FABP4 = "FABP4"
    ACTA2 = "ACTA2"
    CD44 = "CD44"
    IL6 = "IL6"

    # Slide 3
    FOLR1 = "FOLR1"
    CD163 = "CD163"
    TGFB1 = "TGFB1"
    CD8A = "CD8A"

    # Slide 4
    BTLA = "BTLA"
    NCAM1 = "NCAM1"
    IL10 = "IL10"
    MS4A1 = "MS4A1"

    # Slide 5
    NCBP2AS2 = "NCBP2AS2"
    NOTCH3 = "NOTCH3"
    VEGFA = "VEGFA"
    CD34 = "CD34"

    # Slide 6
    IQGAP3 = "IQGAP3"
    CTLA4 = "CTLA4"
    panCK = "panCK"
    CD274 = "CD274"