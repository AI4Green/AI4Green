import numpy as np
import pandas as pd
from scipy.spatial import distance
from sklearn import preprocessing
from sklearn.decomposition import PCA


def getPointChanged(df1, df2):
    df_diff = pd.concat([df1, df2]).drop_duplicates(keep=False)
    try:
        return [
            df_diff.index.tolist()[0],
            [df_diff["PC1"].tolist()[0], df_diff["PC2"].tolist()[0]],
        ]
    except:
        return None
