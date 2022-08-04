import streamlit as st

from CellClasses.RegionType import *
from main_page import PageCreator

page_creator = PageCreator(RegionType.EverythingCells)
st.title("EverythingCells Data")

if "slide_deck" in st.session_state:
    page_creator.create_all()
else:
    st.error("Please load data in the main page first")