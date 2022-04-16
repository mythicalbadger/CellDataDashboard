import streamlit as st

from CellClasses.SlideDeck import *
from PIL import Image

st.set_page_config(layout="wide")


def create_slide_deck(path):
    if "slide_deck" not in st.session_state:
        st.session_state.slide_deck = SlideDeck(path)
    return st.session_state.slide_deck


def get_slide(slide_deck, slide: str):
    slide_num = int(slide.split(" ")[-1])
    return slide_deck.get_slide(slide_num)


def get_tissue(slide, tissue):
    return slide.tissues[tissue]

def save_means(tissue):
    out = tissue.hscores_to_csv()
    return out.to_csv().encode("utf-8")


st.title("Cell Data Dashboard")
st.subheader("Specify Data Path")
st.write("Please specify the **absolute path** to the folder containing your data. The folder hierarchy / naming should look like this:")
st.image(Image.open("resources/folder-hierarchy.png"))
st.write("All data will be in a folder called `Data`. Inside `Data`, there will be a folder called `Slides` for each slide's folder, and a `Thresholds.xlsx` file for the thresholds.")


path = st.text_input(label="Path")
load_btn = st.button(label="Load")

if load_btn:

    deck: SlideDeck = create_slide_deck(path + "Slides/")

    with st.sidebar:
        st.subheader("Choose a slide")
        slide = st.selectbox("Slide", options=["QC Slide", "Slide 1", "Slide 2", "Slide 3", "Slide 4", "Slide 5", "Slide 6"])
        slide = deck.get_slide(slide)

        st.subheader("Choose a tissue")
        tissue = st.selectbox("Tissue", options=[tissue for tissue in slide.tissues])
        tissue = get_tissue(slide, tissue)

        st.subheader("Save H-Score means to CSV")
        csv = save_means(tissue)
        st.download_button(label="Download",
                           data=csv,
                           file_name="hscoremeans.csv",
                           mime="text/csv")

    st.subheader("Cell Counts")
    st.dataframe(tissue.calculate_cell_counts())

    st.subheader("ACD Scores")
    st.dataframe(tissue.score_acd_ranges())

    st.subheader("Zero Scores")
    st.dataframe(tissue.calculate_zero_scores())

    st.subheader("Percent Bins")
    st.dataframe(tissue.score_percent_bins())

    st.subheader("Weighted Percent Bins")
    st.dataframe(tissue.score_weighted_percent_bins())

    st.subheader("H-Scores")
    st.dataframe(tissue.calculate_hscores())

    st.subheader("H-Score Means (Combined)")
    cols = st.columns(4)

    hscore_means = tissue.calculate_hscore_means()
    for i, mean in enumerate(hscore_means):
        cols[i].metric(hscore_means.index[i][0], round(mean, 2))

    hscore_means = tissue.calculate_border_hscore_means()
    if not hscore_means.isnull().values.any():
        st.subheader("H-Score Means (Border)")
        cols = st.columns(4)

        for i, mean in enumerate(hscore_means):
            cols[i].metric(hscore_means.index[i][0], round(mean, 2))

    hscore_means = tissue.calculate_middle_hscore_means()
    if not hscore_means.isnull().values.any():
        st.subheader("H-Score Means (Middle)")
        cols = st.columns(4)

        for i, mean in enumerate(hscore_means):
            cols[i].metric(hscore_means.index[i][0], round(mean, 2))