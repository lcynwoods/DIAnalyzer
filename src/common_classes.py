import logging
import math
import os
import json
import re
import subprocess
import multiprocessing
import statistics
import textwrap
from itertools import combinations
from pathlib import Path

import yaml
import requests
import colorsys
import fastexcel

import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.patches as mpatches
import matplotlib.colorbar as Colorbar
from matplotlib.lines import Line2D

import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.colors as pc
import plotly.io as io
from plotly.subplots import make_subplots

from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import PowerTransformer
from sklearn.metrics import silhouette_score
from sklearn.model_selection import ParameterGrid
from scipy.cluster.hierarchy import linkage, fcluster
import scipy.spatial.distance as ssd
from statsmodels.distributions.empirical_distribution import ECDF
from statsmodels.stats.multitest import multipletests
from yellowbrick.cluster import SilhouetteVisualizer

import networkx as nx
import gseapy as gp
from gseapy import barplot

# Custom modules
from src.Rrvgo_submit import Rrvgo_submit
from src.GO_dot_plot import GO_term_dot_plot
from src.R_normalization_and_analysis import R_normalization_and_analysis
from src.DAPAR_op_intensity_histogram_class import PlotNormCurve
from src.GetDAClass_edit import GetDA
from src.string_networkx import String_network_builder
from src.Kmeans_clustering import KMeansClustering
from src.Single_prot_plotter import POI_plotter
from src.UMAP_cluster import UMAPComparison
from src.PCA_scatter import PCAComparison
from src.hierachical_cluster import HierarchicalCluster
from src.enrichr_api import EnrichrORA
from src.prerank import GSEAPy_Prerank
from src.MA_plotter import PlotMA
from src.Rank_plotter import PlotRank
from src.prepare_data_with_polars import PrepareData
from src.condition_sort import ConditionSort
from src.overview_analysis_with_polars import OverviewAnalysis
from src.prepare_quantification_data import PrepareQuantificationData
from src.GetDAClass_edit import GetDA
from src.push_to_github import upload_html_to_github
