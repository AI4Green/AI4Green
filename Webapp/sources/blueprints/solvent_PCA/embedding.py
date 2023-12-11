import os
from typing import List

import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.decomposition import PCA

from .interactive.Embedder import cPCA
from .interactive.interactive_embeddings import getkPCA


def get_r_class_descriptors(r_class: str) -> (List, List):
    """
    Get descriptors for each reaction class included in Solvent Surfer

    Args:
        r_class: str, reaction class to get descriptors for (case-sensitive)

    Returns:
        class_descs: List of descriptors corresponding to reaction class
        all_descs: List of all descriptors included in dataset
    """
    all_descs = [
        "Molecular Weight",
        "Boiling Point",
        "Density",
        "Viscosity",
        "Vapour Pressure",
        "Refractive Index",
        "LogP",
        "Dipole Moment",
        "Dielectric Constant",
        "Alpha",
        "Beta",
        "Pi",
        "Dispersion",
        "Polarity",
        "H Bonding",
        "Molar Volume",
    ]

    set1 = ["Buchwald-Hartwig coupling"]
    set1_descs = [
        "Dipole Moment",
        "Dielectric Constant",
        "Alpha",
        "Pi",
        "Polarity",
        "H Bonding",
    ]
    set2 = ["Heck cross-coupling", "Grignard", "SNAr/SN2"]
    set2_descs = [
        "Dipole Moment",
        "Dielectric Constant",
        "Alpha",
        "Beta",
        "Pi",
        "Polarity",
        "H Bonding",
    ]
    set3 = ["Alkene metathesis"]
    set3_descs = ["Alpha", "Beta", "Pi", "Dispersion", "H Bonding"]

    if r_class in set1:
        class_descs = set1_descs

    elif r_class in set2:
        class_descs = set2_descs

    elif r_class in set3:
        class_descs = set3_descs

    else:
        class_descs = all_descs

    return class_descs, all_descs


def prepare_data(
    r_class: str,
) -> (
    np.ndarray,
    pd.DataFrame,
    pd.DataFrame,
    List,
    pd.Series,
    pd.Series,
    pd.Series,
    pd.Series,
    pd.Series,
    pd.Series,
    List,
):
    """
    This function prepares the data for kPCA. It extracts the appropriate descriptors for a given reaction class and
    standardises the data using sklearn.StandardScalar(). Returns multiple variables needed for plotting the PCA graphs.

    Args:
        r_class: str, reaction class to generate PCA for

    Returns:
        data_scaled: np.ndarray, scaled data for input to PCA
        class_data_df: pd.DataFrame, contains all descriptors required for reaction class specific PCA
        full_data_df: pd.Dataframe, contains all descriptors from initial dataset
        all_descs: List, contains all descriptors from full_data_df
        names: pd.Series, names of all solvents included in PCA
        CHEM21: pd.Series, CHEM21 classifications of all solvents included in PCA
        CHEM21_numerical: pd.Series, numerical CHEM211 score for each solvent included in PCA
        cost: pd.Series, cost of each solvent included in PCA
        MP: pd.Series, melting point of each solvent included in PCA
        BP: pd.Series, boiling point of each solvent included in PCA
        exp_data: pd.Series, yield data for each reaction class included in PCA
        class_descs: List, reaction specific descriptors used to generate PCA
    """
    # get data from file
    data = pd.read_csv(
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "data", "solvent_data_PCA.csv"
        )
    )
    # descriptors to include in PCA
    class_descs, all_descs = get_r_class_descriptors(r_class)
    class_data_df = data[class_descs]
    full_data_df = data[all_descs]

    # scale data
    scaler = preprocessing.StandardScaler().fit(class_data_df)
    data_scaled = scaler.transform(class_data_df)
    names = data["Solvent"]
    CHEM21 = data["CHEM21"]
    CHEM21_numerical = data["CHEM21 Numerical"]
    cost = data["cost"]
    MP = data["Melting Point"]
    BP = data["Boiling Point"]
    exp_data = data[
        [
            "Grignard",
            "Amide Coupling 1",
            "Amide Coupling 2",
            "Heck",
            "Baylis-Hillman",
            "Suzuki-Miyaura Pd",
            "Suzuki-Miyaura Ni",
            "Buchwald-Hartwig",
            "Alkene Metathesis",
            "SNAr",
        ]
    ]
    return (
        data_scaled,
        class_data_df,
        full_data_df,
        all_descs,
        names,
        CHEM21,
        CHEM21_numerical,
        cost,
        MP,
        BP,
        exp_data,
        class_descs,
    )


def get_PCA(r_class: str) -> (pd.DataFrame, List, pd.Series, cPCA):
    """
    This function gets the initial kernel PCA data for the solvent surfer using a Gaussian kernel.

    Args:
        r_class: str, reaction class to generate kPCA for

    Returns:
        df: pd.DataFrame, dataframe of kPCA data for plotting
        all_descriptors: List, contains all descriptors (not just class specific descriptors)
        names: pd.Series, Contains all name sof solvents included in kPCA
        kPCA: cPCA Object, contains initial embeddings for PCA graph

    """
    (
        data_scaled,
        data_raw,
        data_all,
        all_descriptors,
        names,
        CHEM21,
        CHEM21_numerical,
        cost,
        MP,
        BP,
        exp_data,
        descs,
    ) = prepare_data(r_class)

    # set up kPCA
    PCA_df, kPCA = getkPCA(data_raw, {})

    df = pd.concat(
        [PCA_df, data_all, CHEM21, CHEM21_numerical, names, cost, exp_data], axis=1
    )
    df["MP"] = MP
    df["BP"] = BP
    df = df.fillna("")

    return df, all_descriptors, names, kPCA
