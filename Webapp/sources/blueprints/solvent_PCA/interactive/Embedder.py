#!/usr/bin/python
from collections import defaultdict
from copy import copy

import numpy as np

# from PCA import PCA
from sklearn import decomposition, manifold

# from scipy.spatial import distance
# from PyQt4.QtCore import *
# from PyQt4.QtGui import *
# from PyQt4.Qt import *
# explicitly imported "hidden imports" for pyinstaller
# from sklearn.utils import weight_vector, lgamma
from sklearn.metrics.pairwise import euclidean_distances, pairwise_distances

# Dinos solver
from .cpca import kernel_gen as kernel_gen
from .cpca import solvers as solvers
from .cpca import utils as utils

try:
    from sklearn.neighbors import typedefs
    from sklearn.utils.sparsetools import _graph_validation
except:
    pass


class Embedding(object):
    def __init__(self, data, points, parent):
        self.data = data
        self.original_control_points = None
        self.original_control_point_indices = None
        self.control_points = None
        self.control_point_indices = None
        self.parent = parent
        self.X = np.array([])
        self.Y = np.array([])
        self.ml = []
        self.cl = []
        self.has_ml_cl_constraints = False
        self.projection_matrix = np.zeros((2, len(self.data[0])))
        self.name = ""
        self.is_dynamic = False
        self.update_control_points(points)

    def get_embedding(self):
        pass

    def augment_control_points(self, e):
        avg_median = np.average(abs(np.median(e, axis=0)))
        tmp_points = defaultdict(list)
        if len(self.cl) > 0:
            for pair in self.cl:
                if len(pair) == 2:
                    i, j = list(pair)
                    x1 = e[i]
                    x2 = e[j]
                    diff = x1 - x2
                    norm = np.linalg.norm(diff)
                    new_x1 = x1 + (diff / norm) * 5 * avg_median
                    new_x2 = x2 - (diff / norm) * 5 * avg_median
                    if i not in self.control_point_indices:
                        e[i] = new_x1
                        tmp_points[i] = new_x1
                    if j not in self.control_point_indices:
                        e[j] = new_x2
                        tmp_points[j] = new_x2
        if len(self.ml) > 0:
            for pair in self.ml:
                if len(pair) == 2:
                    i, j = list(pair)
                    x1 = e[i]
                    x2 = e[j]
                    diff = x1 - x2
                    new_x1 = x1 - 0.45 * diff
                    new_x2 = x2 + 0.45 * diff
                    if i not in self.control_point_indices:
                        e[i] = new_x1
                        tmp_points[i] = new_x1
                    if j not in self.control_point_indices:
                        e[j] = new_x2
                        tmp_points[j] = new_x2
        for k, v in tmp_points.items():
            self.control_point_indices.append(k)
            self.control_points.append(v)
        self.X = self.data[self.control_point_indices]
        self.Y = np.array(self.control_points)

    def update_control_points(self, points):
        self.control_point_indices = []
        self.control_points = []
        for i, coords in points.items():
            self.control_point_indices.append(i)
            self.control_points.append(coords)
        self.X = self.data[self.control_point_indices]
        self.Y = np.array(self.control_points)

    def finished_relocating(self):
        pass


class cPCA(Embedding):
    def __init__(self, data, points, parent):
        self.data = data
        self.control_points = []
        self.control_point_indices = []
        self.parent = parent
        self.X = None
        self.Y = np.array([])
        self.projection_matrix = np.zeros((2, len(self.data[0])))
        self.name = ""
        self.is_dynamic = False

        self.ml = []
        self.cl = []
        self.has_ml_cl_constraints = False

        self.name = "cPCA"
        self.projection = np.zeros((2, len(data)))
        self.pca_projection = np.zeros((2, len(data)))
        self.is_dynamic = True
        self.old_control_point_indices = []

        self.params = {
            "r": 3.0,
            "slv_mode": "secular",
            "sigma": None,
            "epsilon": 0.5,
            "degree": 1,
        }
        self.params["const_nu"] = 5e3
        self.params["orth_nu"] = 5e3
        self.params["sigma"] = utils.median_pairwise_distances(data)
        gk = kernel_gen.gaussian_kernel()
        # gk = kernel_gen.polynomial_kernel()
        K = gk.compute_matrix(data, self.params)
        self.embedder = solvers.embedder(2.56e-16, 800, False)
        self.kernel_sys = self.embedder.kernel_sys(K)
        # self.parent.status_text.setText("Done, calculating Gaussean kernel.")

        label_mask = np.array([0])
        self.quad_eig_sys = self.embedder.sph_cl_var_term_eig_sys(self.kernel_sys)
        self.quad_eig_sys_original = copy(self.quad_eig_sys)
        if len(self.control_point_indices) == 0:
            placement_mask = np.array([0])
        else:
            placement_mask = np.array(self.control_point_indices)
        self.const_mu = self.embedder.const_nu(
            self.params, placement_mask, self.kernel_sys
        )
        self.update_control_points(points)
        self.finished_relocating()
        if len(self.Y) == 0:
            pca_dirs = self.embedder.soft_cp_mode_directions(
                self.quad_eig_sys,
                label_mask,
                np.ones((1, 2)),
                self.kernel_sys,
                self.params,
                1e-20,
            )
        else:
            for i in range(len(self.control_point_indices)):
                self.quad_eig_sys = self.embedder.sph_cp_quad_term_eig_sys(
                    self.kernel_sys,
                    self.quad_eig_sys,
                    self.control_point_indices[i],
                    self.const_mu,
                )
            pca_dirs = self.embedder.soft_cp_mode_directions(
                self.quad_eig_sys,
                self.control_point_indices,
                self.Y,
                self.kernel_sys,
                self.params,
                self.const_mu,
            )
        self.pca_projection = self.kernel_sys[0].dot(pca_dirs)

    def get_embedding(self, X=None):
        if set(self.control_point_indices) != self.old_control_point_indices:
            self.pca_projection = self.finished_relocating()
        self.old_control_point_indices = set(self.control_point_indices)
        return self.pca_projection.T

    def finished_relocating(self):
        if len(self.control_point_indices) > 0:
            directions = self.embedder.soft_cp_mode_directions(
                self.quad_eig_sys,
                self.control_point_indices,
                self.Y,
                self.kernel_sys,
                self.params,
                self.const_mu,
            )
            self.pca_projection = self.kernel_sys[0].dot(directions)
        return self.pca_projection

    def update_control_points(self, points):
        super(cPCA, self).update_control_points(points)
        if len(self.control_point_indices) > len(self.old_control_point_indices):
            selected_point = None
            if selected_point is None:
                selected_point = (
                    list(
                        set(self.control_point_indices)
                        - set(self.old_control_point_indices)
                    )
                )[0]

            self.quad_eig_sys = self.embedder.sph_cp_quad_term_eig_sys(
                self.kernel_sys, self.quad_eig_sys, selected_point, self.const_mu
            )
            directions = self.embedder.soft_cp_mode_directions(
                self.quad_eig_sys,
                self.control_point_indices,
                self.Y,
                self.kernel_sys,
                self.params,
                self.const_mu,
            )
            self.pca_projection = self.kernel_sys[0].dot(directions)
        elif len(self.control_point_indices) < len(self.old_control_point_indices):
            self.quad_eig_sys = copy(self.quad_eig_sys_original)
            for i in range(len(self.control_point_indices)):
                self.quad_eig_sys = self.embedder.sph_cp_quad_term_eig_sys(
                    self.kernel_sys,
                    self.quad_eig_sys,
                    self.control_point_indices[i],
                    self.const_mu,
                )
            if len(self.control_point_indices) == 0:
                pca_dirs = self.embedder.soft_cp_mode_directions(
                    self.quad_eig_sys,
                    np.array([0]),
                    np.ones((1, 2)),
                    self.kernel_sys,
                    self.params,
                    1e-20,
                )
                self.pca_projection = self.kernel_sys[0].dot(pca_dirs)
            else:
                pca_dirs = self.embedder.soft_cp_mode_directions(
                    self.quad_eig_sys,
                    self.control_point_indices,
                    self.Y,
                    self.kernel_sys,
                    self.params,
                    self.const_mu,
                )
                self.pca_projection = self.kernel_sys[0].dot(pca_dirs)
        self.old_control_point_indices = set(self.control_point_indices)

        if self.has_ml_cl_constraints:
            self.augment_control_points(self.get_embedding().T)
            self.quad_eig_sys = copy(self.quad_eig_sys_original)
            for i in range(len(self.control_point_indices)):
                self.quad_eig_sys = self.embedder.sph_cp_quad_term_eig_sys(
                    self.kernel_sys,
                    self.quad_eig_sys,
                    self.control_point_indices[i],
                    self.const_mu,
                )
            if len(self.control_point_indices) == 0:
                pca_dirs = self.embedder.soft_cp_mode_directions(
                    self.quad_eig_sys,
                    np.array([0]),
                    np.ones((1, 2)),
                    self.kernel_sys,
                    self.params,
                    1e-20,
                )
                self.pca_projection = self.kernel_sys[0].dot(pca_dirs)
            else:
                pca_dirs = self.embedder.soft_cp_mode_directions(
                    self.quad_eig_sys,
                    self.control_point_indices,
                    self.Y,
                    self.kernel_sys,
                    self.params,
                    self.const_mu,
                )
                self.pca_projection = self.kernel_sys[0].dot(pca_dirs)
