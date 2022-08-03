import streamlit as st

from CellClasses.SlideDeck import *
from CellClasses.RegionType import *
from main_page import PageCreator

page_creator = PageCreator(RegionType.EverythingCells)
st.title("EverythingCells Data")
def get_tissue(slide, tissue):
    return slide.tissues[tissue]

if "slide_deck" in st.session_state:
    page_creator.create_all()
else:
    st.markdown("Please load data in the main page first")