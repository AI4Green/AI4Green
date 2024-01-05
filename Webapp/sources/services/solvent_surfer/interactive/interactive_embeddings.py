import numpy as np
import pandas as pd

from .Dataset import Dataset
from .Embedder import cPCA


class GetEmbedding:
    def __init__(self, pandasdata=""):
        super(GetEmbedding, self).__init__()
        self.embedding = None
        self.embedding_algorithm = None
        self.data = Dataset()
        self.data.load_pandas_dataframe(pandasdata)

    def update(self):
        """Renders the embedding"""
        self.embedding = self.embedding_algorithm.get_embedding()
        return np.array([self.embedding[0], self.embedding[1]])

    def load_embedding(self):
        """
        renders embedding from preloaded embedding algorithm
        """
        self.embedding = self.embedding_algorithm.pca_projection.T
        return np.array([self.embedding[0], self.embedding[1]])

def getkPCA(data, mode, embedding=None):
    kPCA = GetEmbedding(pandasdata=data)

    kPCA.embedding_algorithm = cPCA(kPCA.data.data, {}, kPCA)

    if mode == 'on_load':

        kPCA.embedding_algorithm.load_attrs(embedding)
        points = kPCA.load_embedding()

    else:
        points = kPCA.update()

    df = pd.DataFrame(data=points.T, columns=["PC1", "PC2"])

    return df, kPCA
