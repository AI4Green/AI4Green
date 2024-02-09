"""
Created on May 14, 2013

@author: doglic
"""
import sys

import numpy as np
import scipy.spatial.distance as dist

sys.path.append("./../common")
from . import utils


def center_kernel(K):
    n = K.shape[0]
    H = np.identity(n) - float(1.0 / n) * np.ones(K.shape)
    return H.dot(K).dot(H)


def cos_transform(K):
    v = np.zeros(K.shape[0])
    for i in range(v.shape[0]):
        v[i] = float(1.0 / np.sqrt(K[i, i]))

    return K * np.dot(v.reshape(-1, 1), v.reshape(1, -1))


class kernel_functional(object):
    def compute_matrix(self, points, params=None):
        pass

    def evaluate(self, x, y, params=None):
        pass


class polynomial_kernel(kernel_functional):
    # expected to be n x dim
    def compute_matrix(self, points, params={"degree": 1}):
        M = np.dot(points, np.transpose(points))
        return np.power(1 + M, params["degree"])

    def evaluate(self, x, y, params={"degree": 1}):
        x = x.reshape(1, -1)
        y = y.reshape(-1, 1)

        return np.power(1 + np.dot(x, y), params["degree"])


class cos_linear_kernel(polynomial_kernel):
    def compute_matrix(self, points, params={"degree": 1}):
        M = np.dot(points, np.transpose(points))
        M = np.power(1 + M, params["degree"])
        return cos_transform(M)

    def evaluate(self, x, y, params={"degree": 1}):
        raise ValueError


class gaussian_kernel(kernel_functional):
    def compute_matrix(self, points, params={"sigma": 1.0}):
        euclidean_distances = dist.squareform(dist.pdist(points, "euclidean"))
        gauss_distances = np.exp(
            -(euclidean_distances**2) / (2 * params["sigma"] ** 2)
        )
        return gauss_distances


class pinv_k_nn_laplacian_kernel(kernel_functional):
    def compute_matrix(self, points, params={"sigma": 1.0, "k": 3}):
        adj_mat = float(
            0.5 / np.power(params["sigma"], 2)
        ) * utils.k_nn_pw_sq_distance_mat(points, params["k"])
        adj_mat = np.exp(-adj_mat)

        adj_mat_normalizer = np.diag(1.0 / np.sqrt(np.sum(adj_mat, axis=0)))
        laplacian_mat = np.identity(adj_mat.shape[0]) - np.dot(
            np.dot(adj_mat_normalizer, adj_mat), adj_mat_normalizer
        )

        pinv_laplacian_mat = np.linalg.pinv(laplacian_mat) + 1e-6 * np.identity(
            adj_mat.shape[0]
        )

        # make sure it is symmetric and tweak numerical errors
        return 0.5 * (pinv_laplacian_mat + np.transpose(pinv_laplacian_mat))

    def evaluate(self, x, y, params=None):
        raise ValueError(
            "Not supported! You must compute the full kernel matrix as the dot product depends on all data points!"
        )


class isomap_kernel(kernel_functional):
    def compute_matrix(self, points, params=None):
        n = points.shape[0]
        # step I: compute the nn graph
        edist_adj_mat = utils.pairwise_distance_mat(points)
        g = utils.sim_mat_to_graph(edist_adj_mat)
        mst_result = utils.min_spanning_tree(g, 0)
        cutoff = utils.largest_mst_weigth(points, mst_result[0])
        g = utils.sim_mat_to_graph(edist_adj_mat, cutoff, False)

        # step II: compute the pairwise distances using Dijkstra algorithm based on the computed adjacency matrix
        shortest_path_mat = utils.adj_mat_to_dijkstra_spath_mat(g)

        # step III: convert to the Gram matrix
        H = np.identity(n) - float(1.0 / n) * np.ones(shortest_path_mat.shape)
        return -0.5 * np.dot(np.dot(H, shortest_path_mat), H)

    def evaluate(self, x, y, params=None):
        raise ValueError(
            "Not supported! You must compute the full kernel matrix as the dot product depends on all data points!"
        )


class pinv_laplacian_kernel(kernel_functional):
    def compute_matrix(self, points, params=None):
        n = points.shape[0]
        edist_adj_mat = utils.pairwise_distance_mat(points)
        if params is None or params["sigma"] is None:
            g = utils.sim_mat_to_graph(edist_adj_mat)
            mst_result = utils.min_spanning_tree(g, 0)
            gauss_sigma = utils.largest_mst_weigth(points, mst_result[0])
        else:
            gauss_sigma = params["sigma"]

        gauss_adj_mat = utils.prune_kernel_mat(
            np.exp(float(-0.5 / np.power(gauss_sigma, 2)) * np.power(edist_adj_mat, 2)),
            params["epsilon"],
        )
        adj_deg_normalizer = np.diag(1.0 / np.sqrt(np.sum(gauss_adj_mat, axis=0)))
        laplacian_mat = np.identity(n) - np.dot(
            np.dot(adj_deg_normalizer, gauss_adj_mat), adj_deg_normalizer
        )
        pinv_laplacian = np.linalg.pinv(laplacian_mat) + 2.56e-6 * np.eye(n)

        # make sure it is symmetric and tweak numerical errors
        return 0.5 * (pinv_laplacian + np.transpose(pinv_laplacian))

    def evaluate(self, x, y, params=None):
        raise ValueError(
            "Not supported! You must compute the full kernel matrix as the dot product depends on all data points!"
        )
