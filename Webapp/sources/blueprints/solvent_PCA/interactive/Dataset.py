#!/usr/bin/python
import os
import sys
from copy import copy

import numpy as np
import pandas as pd
from sklearn.datasets import load_svmlight_file

try:
    from sklearn.utils import (  # needed as explicit import for generating one binary executable file
        arraybuilder,
    )
except:
    pass
from statsmodels.tsa.ar_model import AR
from statsmodels.tsa.arima_model import ARMA


class Dataset:
    def __init__(self):
        self.database_type = None
        self.original_data = []
        self.data = []
        self.instance_names = []
        self.attribute_names = []
        self.ignored_attributes = []
        self.considered_attributes = []
        self.label_name = ""
        self.label_index = -1
        self.dataset_name = ""
        self.masked_names = None
        self.pictures = []
        self.raw_data = []
        self.normalize = self.normalize_max

    def load_pandas_dataframe(self, data):
        self.dataset_name = "Pandas DataFrame"
        self.database_type = "DATA"
        self.attribute_names = np.array([str(i) for i in data.columns])
        if data.index.dtype.kind in "biufc":  # is_numeric()
            self.instance_names = ["Instance: " + str(s) for s in data.index]
        else:
            self.instance_names = list(data.index.astype(str))
        self.considered_attributes = (
            np.ones(len(self.attribute_names)).astype(int).tolist()
        )
        self.original_data = np.array(data.values).astype(float)

        stds = np.std(self.original_data, axis=0)
        keeper = list(np.nonzero(stds)[0])
        self.attribute_names = self.attribute_names[keeper]
        self.original_data = self.original_data.T[keeper].T

        self.data = copy(self.original_data)
        self.label_name = self.attribute_names[self.label_index]
        self.ignored_attributes.append(self.label_name)
        self.mask = []
        for i, name in enumerate(self.attribute_names):
            if name not in self.ignored_attributes:
                self.mask.append(i)
        self.update_data()
        self.normalize()

    def normalize_max(self, X=[]):
        if X == []:
            self.mask = []
            for i, name in enumerate(self.attribute_names):
                if name not in self.ignored_attributes:
                    self.mask.append(i)
            self.masked_names = [self.attribute_names[i] for i in self.mask]
            self.data = copy(self.original_data.T[self.mask].T)
            means = np.mean(self.data, axis=0)
            shifted_data = self.data - means
            maxval = np.max(shifted_data, axis=0)
            self.data = shifted_data / maxval
        else:
            means = np.mean(X, axis=0)
            shifted_data = X - means
            maxval = np.max(shifted_data, axis=0)
            X = shifted_data / maxval
            return X

    def update_data(self):
        self.mask = []
        for i, name in enumerate(self.attribute_names):
            if name not in self.ignored_attributes:
                self.mask.append(i)
        self.masked_names = [self.attribute_names[i] for i in self.mask]
        self.data = copy(self.original_data.T[self.mask].T)
        self.normalize()
        if np.isnan(self.data).any():
            print("Nan in copied data after normalizing:")
            print(self.original_data)
