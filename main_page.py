import streamlit as st

from CellClasses.SlideDeck import *
from CellClasses.RegionType import *
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


class PageCreator():
    def __init__(self, region: RegionType):
        self.region = region
        self.tissue = None

    def create_download_button(self, name: str, data: pd.DataFrame):
        st.download_button(label="Download",
                           data=data.to_csv(encoding="utf-8"),
                           file_name=self.tissue.name + f"-{name}-{self.region.value}.csv",
                           mime="text/csv",
                           key=f"{self.tissue.name}_{name}_{self.region.value}")

    def create_hscore_means(self):
        return self.tissue.hscores_to_csv(self.region).to_csv().encode("utf-8")

    def create_all(self):
        deck: SlideDeck = st.session_state.slide_deck
        with st.sidebar:
            st.subheader("Choose a slide")
            slide = st.selectbox("Slide",
                                 options=["QC Slide", "Slide 1", "Slide 2", "Slide 3", "Slide 4", "Slide 5", "Slide 6"],
                                 key=f"selectbox-{self.region.value}")
            slide = deck.get_slide(slide)

            st.subheader("Choose a tissue")
            tissue = st.selectbox("Tissue", options=[tissue for tissue in slide.tissues],
                                  key=f"selectbox-{self.region.value}")
            tissue = get_tissue(slide, tissue)
            self.tissue = tissue

        st.header("Cell Counts")
        df = tissue.calculate_cell_counts(self.region)
        self.create_download_button("cell-counts", df)
        st.dataframe(df)

        st.header("Cell Expressions")
        df = tissue.calculate_expressed_markers(self.region)
        self.create_download_button("expressed-markers", df)
        st.dataframe(df)

        st.header("ACD Scores")
        df = tissue.score_acd_ranges(self.region)
        self.create_download_button("acd-scores", df)
        st.dataframe(df)

        st.header("Median Markers")
        df = tissue.calculate_median_markers(self.region)
        self.create_download_button("median-markers", df)
        st.dataframe(df)

        st.header("Zero Scores")
        df = tissue.calculate_zero_scores(self.region)
        self.create_download_button("zero-scores", df)
        st.dataframe(df)

        st.header("Positive Expressors")
        df = tissue.calculate_positive_expressors(self.region)
        self.create_download_button("positive-expressors", df)
        st.dataframe(df)

        st.header("Percent Bins")
        df = tissue.score_percent_bins(self.region)
        self.create_download_button("percent-bins", df)
        st.dataframe(df)

        st.header("Weighted Percent Bins")
        df = tissue.score_weighted_percent_bins(self.region)
        self.create_download_button("weighted-percent-bins", df)
        st.dataframe(df)

        st.header("Median Neighbors")
        df = tissue.median_nbrs(self.region)
        self.create_download_button("median-neighbors", df)
        st.dataframe(df)

        st.header("Median Percent Touching")
        df = tissue.median_percent_touching(self.region)
        self.create_download_button("median-percent-touching", df)
        st.dataframe(df)

        st.header("Pair Gene Expression")
        df = tissue.calculate_pair_gene_expression(self.region)
        self.create_download_button("pair-gene-expression-touching", df)
        st.dataframe(df)

        st.header("Trio Gene Expression")
        df = tissue.calculate_trio_gene_expression(self.region)
        self.create_download_button("trio-gene-expression-touching", df)
        st.dataframe(df)

        st.header("Quad Gene Expression")
        df = tissue.calculate_quad_gene_expression(self.region)
        self.create_download_button("quad-gene-expression-touching", df)
        st.dataframe(df)

        st.header("H-Scores")
        df = tissue.calculate_hscores(self.region)
        self.create_download_button("hscores", df)
        st.dataframe(df)

        st.header("H-Score Medians (Combined)")
        st.download_button(label="Download",
                           data=self.create_hscore_means(),
                           file_name=tissue.name + "-hscoremeans.csv",
                           mime="text/csv")
        cols = st.columns(4)

        hscore_means = tissue.calculate_hscore_medians(self.region)
        for i, mean in enumerate(hscore_means):
            cols[i].metric(hscore_means.index[i][0], round(mean, 2))

        hscore_means = tissue.calculate_border_hscore_medians(self.region)
        if not hscore_means.isnull().values.any():
            st.header("H-Score Medians (Border)")
            cols = st.columns(4)

            for i, mean in enumerate(hscore_means):
                cols[i].metric(hscore_means.index[i][0], round(mean, 2))

        hscore_means = tissue.calculate_middle_hscore_medians(self.region)
        if not hscore_means.isnull().values.any():
            st.header("H-Score Medians (Middle)")
            cols = st.columns(4)

            for i, mean in enumerate(hscore_means):
                cols[i].metric(hscore_means.index[i][0], round(mean, 2))


st.title("Cell Data Dashboard")
st.subheader("Specify Data Path")
st.write(
    "Please specify the **absolute path** to the folder containing your data. The folder hierarchy / naming should look like this:")
st.image(Image.open("resources/folder-hierarchy.png"))
st.write(
    "All data will be in a folder called `Data`. Inside `Data`, there will be a folder called `Slides` for each slide's folder, and a `Thresholds.xlsx` file for the thresholds.")

path = st.text_input(label="Path", value="/home/mythicalbadger/CellDataExtractor/Data/")
load_btn = st.button(label="Load")

if load_btn and "slide_deck" not in st.session_state:
    deck: SlideDeck = create_slide_deck(path)

if "slide_deck" in st.session_state:
    deck: SlideDeck = st.session_state.slide_deck
    with st.sidebar:
        st.subheader("Choose a slide")
        slide = st.selectbox("Slide",
                             options=["QC Slide", "Slide 1", "Slide 2", "Slide 3", "Slide 4", "Slide 5", "Slide 6"])
        slide = deck.get_slide(slide)

        st.subheader("Choose a tissue")
        tissue = st.selectbox("Tissue", options=[tissue for tissue in slide.tissues])
        tissue = get_tissue(slide, tissue)
    st.success("Data loaded successfully")
