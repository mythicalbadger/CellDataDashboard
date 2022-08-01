from CellClasses.SlideDeck import *

deck: SlideDeck = SlideDeck("/home/mythicalbadger/CellDataExtractor/Data/")
slide = deck.get_slide("Slide 1")
tissue = slide.tissues["ZCS008 series"]
print(tissue.median_nbrs())