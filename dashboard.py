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

def save_hscores(tissue):
    out = tissue.means_to_csv()
    return out.to_csv().encode("utf-8")

def save_zero_scores(tissue):
    return tissue.calculate_zero_scores().to_csv().encode("utf-8")


st.title("Cell Data Dashboard")
st.subheader("Specify Data Path")
st.write("Please specify the **absolute path** to the folder containing your data. The folder hierarchy / naming should look like this:")
st.image(Image.open("resources/folder-hierarchy.png"))
st.write("All data will be in a folder called `Data`. Inside `Data`, there will be a folder called `Slides` for each slide's folder, and a `Thresholds.xlsx` file for the thresholds.")


path = st.text_input(label="Path")
load_btn = st.button(label="Load")

if load_btn and "slide_deck" not in st.session_state:
    deck: SlideDeck = create_slide_deck(path)

if "slide_deck" in st.session_state:
    deck: SlideDeck = st.session_state.slide_deck
    with st.sidebar:
        st.subheader("Choose a slide")
        slide = st.selectbox("Slide", options=["QC Slide", "Slide 1", "Slide 2", "Slide 3", "Slide 4", "Slide 5", "Slide 6"])
        slide = deck.get_slide(slide)

        st.subheader("Choose a tissue")
        tissue = st.selectbox("Tissue", options=[tissue for tissue in slide.tissues])
        tissue = get_tissue(slide, tissue)

    st.subheader("Cell Counts")
    cell_counts = tissue.calculate_cell_counts()
    st.download_button(label="Download",
                       data=cell_counts.to_csv().encode("utf-8"),
                       file_name=tissue.name + "-cell-counts.csv",
                       mime="text/csv")
    st.dataframe(cell_counts)

    st.subheader("Cell Expressions")
    expressed_markers = tissue.calculate_expressed_markers()
    st.download_button(label="Download",
                       data=expressed_markers.to_csv().encode("utf-8"),
                       file_name=tissue.name + "-expressed-markers.csv",
                       mime="text/csv")
    st.dataframe(expressed_markers)

    st.subheader("Median Markers")
    median_markers = tissue.calculate_median_markers()
    st.download_button(label="Download",
                       data=median_markers.to_csv().encode("utf-8"),
                       file_name=tissue.name + "-median-markers.csv",
                       mime="text/csv")
    st.dataframe(median_markers)

    st.subheader("Zero Scores")
    zero_scores = tissue.calculate_zero_scores()
    st.download_button(label="Download",
                       data=zero_scores.to_csv().encode("utf-8"),
                       file_name=tissue.name + "-zero_scores.csv",
                       mime="text/csv")
    st.dataframe(zero_scores)

    st.subheader("ACD Scores")
    acd_ranges = tissue.score_acd_ranges()
    st.download_button(label="Download",
                       data=acd_ranges.to_csv().encode("utf-8"),
                       file_name=tissue.name + "-acd_ranges.csv",
                       mime="text/csv")
    st.dataframe(acd_ranges)

    st.subheader("Positive Expressors")
    positive_expressors = tissue.calculate_positive_expressors()
    st.download_button(label="Download",
                       data=positive_expressors.to_csv().encode("utf-8"),
                       file_name=tissue.name + "-positive_expressors.csv",
                       mime="text/csv")
    st.dataframe(positive_expressors)

    st.subheader("Percent Bins")
    percent_bins = tissue.score_percent_bins()
    st.download_button(label="Download",
                       data=percent_bins.to_csv().encode("utf-8"),
                       file_name=tissue.name + "-percent_bins.csv",
                       mime="text/csv")
    st.dataframe(percent_bins)

    st.subheader("Weighted Percent Bins")
    weighted_percent_bins = tissue.score_weighted_percent_bins()
    st.download_button(label="Download",
                       data=weighted_percent_bins.to_csv().encode("utf-8"),
                       file_name=tissue.name + "-weighted_percent_bins.csv",
                       mime="text/csv")
    st.dataframe(weighted_percent_bins)

    st.subheader("H-Scores")
    hscores = tissue.calculate_hscores()
    st.download_button(label="Download",
                       data=hscores.to_csv().encode("utf-8"),
                       file_name=tissue.name + "-hscores.csv",
                       mime="text/csv")
    st.dataframe(hscores)

    st.subheader("H-Score Means (Combined)")
    st.download_button(label="Download",
                       data=save_means(tissue),
                       file_name=tissue.name + "-hscoremeans.csv",
                       mime="text/csv")
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

    st.subheader("Median Neighbors")
    median_nbrs = tissue.median_nbrs()
    st.download_button(label="Download",
                       data=median_nbrs.to_csv().encode("utf-8"),
                       file_name=tissue.name + "-median-neighbors.csv",
                       mime="text/csv")
    st.dataframe(median_nbrs)
