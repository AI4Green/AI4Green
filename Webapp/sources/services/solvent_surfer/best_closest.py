from math import dist

import numpy as np
import pandas as pd


def get_closest(index: int, data: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates Euclidean distance between point with index and all other data points

    Args:
        index - index corresponding to datapoint selected from PCA graph
        data - kPCA data including data for first two principal components

    Returns:
        pd.DataFrame with Euclidean distances from point with index, organised closest to farthest
    """
    selected_point = data.iloc[index]
    distances = []
    for i in range(len(data["PC1"])):
        distances.append(
            dist(
                (selected_point["PC1"], selected_point["PC2"]),
                (data["PC1"].tolist()[i], data["PC2"].tolist()[i]),
            )
        )
    distances = np.array(distances)
    normalized_distances = (distances - np.min(distances)) / (
        np.max(distances) - np.min(distances)
    )
    data["distances"] = normalized_distances
    data = data.sort_values(by="distances")
    data = data.reset_index(drop=True)

    return data
