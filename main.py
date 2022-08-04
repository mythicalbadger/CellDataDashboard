from CellClasses.SlideDeck import *
from CellClasses.RegionType import *

deck: SlideDeck = SlideDeck("/home/mythicalbadger/CellDataExtractor/Data/")
slide = deck.get_slide("Slide 1")
tissue = slide.tissues["ZCS008 series"]
print(tissue.calculate_pair_gene_expression(RegionType.EverythingCells))