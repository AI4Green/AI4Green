"""
Created on Nov 28, 2013

@author: doglic
"""

import heapq
import logging
import multiprocessing as mproc

# import matplotlib
import warnings
from logging import handlers

# import time
# import os
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial.distance as dist

# import pickle


# def async_min_spanning_tree(*args):
#    mm_g = None
#    while mm_g == None:
#        try:
#            f = open('tmp.dat', 'rb')
#            mm_g = pickle.load(f)
#            f.close()
#        except:
#            time.sleep(1)
#    return min_spanning_tree(mm_g, args[0])


def get_log_file_handler(fileName, maxCapacity, backups, level):
    handler = handlers.RotatingFileHandler(
        fileName, maxBytes=maxCapacity, backupCount=backups
    )
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    handler.setLevel(level)
    return handler


def prune_kernel_mat(K, threshold):
    n = K.shape[0]
    for i in range(n):
        for j in range(i):
            if K[i, j] < threshold:
                K[i, j] = 0.0
                K[j, i] = 0.0

    return K


def sqrt_inverse_sym_pd_mat(K):
    eigvals, eigvecs = np.linalg.eigh(K)

    min_real = float(min(np.real(eigvals)))
    max_abs_imag = float(max(np.abs(np.imag(eigvals))))
    if min_real <= 0 or max_abs_imag > 0:
        warnings.warn(
            "WARNING: The kernel matrix K is poorly conditioned! Adding noise to the diagonal..."
        )
        eigvals = eigvals + np.abs(min_real) + 1.0e-6 * np.random.rand()

    D = np.zeros((eigvals.shape[0], eigvals.shape[0]))
    for i in range(D.shape[0]):
        if eigvals[i] <= 0 or np.abs(np.imag(eigvals[i])) > 0:
            raise ValueError(
                "The matrix K is not strictly positive definite, Eigenvalue_"
                + str(i)
                + " = "
                + str(eigvals[i])
            )

        D[i][i] = float(1.0 / np.sqrt(np.real(eigvals[i])))

    U = np.asmatrix(eigvecs)

    return np.dot(np.dot(U, D), np.transpose(U))


def rect_mat_qr(L):
    m = L.shape[0]
    n = L.shape[1]
    M = np.zeros((n, n))
    M[:m, :] = L[:, :]
    Q, R = np.linalg.qr(np.transpose(M))
    return (Q, R)


def is_soft_orthogonal(logger, alphas, K, threshold=0.95):
    d = alphas.shape[1]
    angles = []
    for i in range(d):
        for j in range(i + 1, d):
            orth_coeff = np.dot(np.dot(alphas[:, i], K), alphas[:, j])
            angles.append(orth_coeff)
            logger.debug(
                "< f_{" + str(i) + "}, f_{" + str(j) + "} > = " + str(orth_coeff)
            )

    if np.mean(np.abs(angles)) < threshold:
        return True

    return False


def is_strictly_soft_orthogonal(logger, alphas, K, threshold=0.95):
    d = alphas.shape[1]
    for i in range(d):
        for j in range(i + 1, d):
            orth_coeff = np.dot(np.dot(alphas[:, i], K), alphas[:, j])
            logger.debug(
                "< f_{" + str(i) + "}, f_{" + str(j) + "} > = " + str(orth_coeff)
            )
            if np.abs(orth_coeff) > threshold:
                return False

    return True


def nn_cl_error(X, y, labelled_idxs, boundary_plot_path=None):
    pred_labels = nn_labels(X, y, labelled_idxs)
    if boundary_plot_path is not None:
        planar_class_embedding(X, pred_labels, labelled_idxs, boundary_plot_path)
    num_errors = np.count_nonzero(pred_labels - y)
    return float(num_errors * 1.0 / max(1.0, y.shape[0] - labelled_idxs.shape[0]))


def nn_labels(X, y, labelled_idxs):
    n = X.shape[0]
    pred_labels = np.zeros((n, 1))
    for i in range(n):
        d = np.inf
        for j in labelled_idxs:
            if i == j:
                pred_labels[i] = y[j]
                break
            else:
                d2j = np.linalg.norm(X[i, :] - X[j, :])
                if d2j < d:
                    pred_labels[i] = y[j]
                    d = d2j

    return pred_labels


def center_data(X):
    mu = np.mean(X, axis=0)
    sigma = np.std(X, axis=0)

    return (X - mu) / sigma


def nn_distances(X):
    n = X.shape[0]
    nn_dists = 1000000 * np.ones(n)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue

            d = np.linalg.norm(X[i, :] - X[j, :])
            if d < nn_dists[i]:
                nn_dists[i] = d

    return nn_dists


def pairwise_distance_mat(X):
    n = X.shape[0]
    pd_adj_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i):
            pd_adj_mat[i, j] = np.linalg.norm(X[i, :] - X[j, :])
            pd_adj_mat[j, i] = pd_adj_mat[i, j]

    return pd_adj_mat


def k_nn_pw_sq_distance_mat(X, k):
    n = X.shape[0]
    adj_mat = np.zeros((n, n))
    for i in range(n):
        top_k_nn_dists = np.Inf * np.ones(k)
        top_k_nn_idxs = np.zeros(k)
        threshold_tuple = (0, np.Inf)
        for j in range(n):
            if i == j:
                continue
            i2j_vec = X[i, :] - X[j, :]
            i2j_dist = float(np.dot(i2j_vec.reshape(1, -1), i2j_vec.reshape(-1, 1)))
            if i2j_dist < threshold_tuple[1]:
                insert_pos = np.argmax(top_k_nn_dists)
                top_k_nn_dists[insert_pos] = i2j_dist
                top_k_nn_idxs[insert_pos] = j

                new_th_idx = np.argmax(top_k_nn_dists)
                threshold_tuple = (new_th_idx, top_k_nn_dists[new_th_idx])

        for j in range(k):
            adj_mat[i, top_k_nn_idxs[j]] = top_k_nn_dists[j]
            adj_mat[top_k_nn_idxs[j], i] = adj_mat[i, top_k_nn_idxs[j]]

    return adj_mat


# def min_spanning_tree(adj_mat, root=0):
#    n = adj_mat.shape[0]
#    c_inf = np.sum(adj_mat)
#    vertex_costs = c_inf * np.ones(n)
#    distances = np.zeros(n)
#    incid_vec = -np.ones(n)
#    visited = np.zeros(n)
#    visited[root] = 1
#    min_cost_vertex = root
#
#    while np.sum(visited) != n:
#        for i in range(n):
#            if visited[i] == 0 and vertex_costs[i] > adj_mat[i, min_cost_vertex]:
#                    vertex_costs[i] = adj_mat[i, min_cost_vertex]
#                    incid_vec[i] = min_cost_vertex
#                    distances[i] = distances[min_cost_vertex] + adj_mat[i, min_cost_vertex]
#
#        min_cost_vertex = np.argmin(vertex_costs)
#        visited[min_cost_vertex] = 1
#        vertex_costs[min_cost_vertex] = c_inf
#
#    return (incid_vec, distances)


"""def sim_mat_to_graph(sim_mat, threshold=0.0, prune_lower=True):
    n = sim_mat.shape[0]
    g = graph.graph(n)
    for i in range(n):
        for j in range(i):
            if prune_lower and sim_mat[i, j] > threshold:
                g.add_edge(i, j, sim_mat[i, j])
            elif prune_lower == False and sim_mat[i, j] < threshold:
                g.add_edge(i, j, sim_mat[i, j])
    return g"""


def min_spanning_tree(g, root=0):
    try:
        n = g.v()
        mst_edges = -np.ones(n)
        sh_paths = np.zeros(n)
        visited = np.zeros(n)
        v_costs = np.Inf * np.ones(n)
        mc_vertex = root
        visited[mc_vertex] = 1

        while True:
            for e in g.adj(mc_vertex):
                if visited[e.head()] == 0 and v_costs[e.head()] > e.weight():
                    v_costs[e.head()] = e.weight()
                    mst_edges[e.head()] = e.tail()
                    sh_paths[e.head()] = sh_paths[e.tail()] + e.weight()

            mc_vertex = np.argmin(v_costs)
            if visited[mc_vertex] == 1:
                break

            visited[mc_vertex] = 1
            v_costs[mc_vertex] = np.Inf

        return (mst_edges, sh_paths)
    except Exception as ex:
        warnings.warn("error while computing minimum spanning tree: " + str(ex))
        return None


# def pq_min_spanning_tree(adj_mat, root):
#    n = adj_mat.shape[0]
#    edge_pq = []
#    distances = np.zeros(n)
#    incid_vec = -np.ones(n)
#    visited = np.zeros(n)
#    return __pq_prim(adj_mat, root, edge_pq, visited, distances, incid_vec)


def pq_min_spanning_tree(g, root=0):
    n = g.v()
    epq = []
    mst_edges = -np.ones(n)
    sh_paths = np.zeros(n)
    visited = np.zeros(n)
    return pq_prim(g, root, epq, visited, mst_edges, sh_paths)


def prim_scan(g, u, epq, visited):
    visited[u] = 1
    for e in g.adj(u):
        if visited[e.head()] == 0:
            heapq.heappush(epq, e)


# def __prim_scan(edge_pq, adj_mat, u, visited):
#    visited[u] = 1
#    for i in range(adj_mat.shape[0]):
#        if adj_mat[u, i] > 0 and visited[i] == 0:
#            heapq.heappush(edge_pq, graph.edge(u, i, adj_mat[u, i]))


def pq_prim(g, root, epq, visited, mst_edges, sh_paths):
    try:
        prim_scan(g, root, epq, visited)
        num_visited = 1
        while num_visited < g.v():
            e = heapq.heappop(epq)
            u = e.tail()
            v = e.head()

            if visited[v] == 0:
                mst_edges[v] = u
                sh_paths[v] = sh_paths[u] + e.weight()
                prim_scan(g, v, epq, visited)
                num_visited += 1
    except IndexError:
        pass

    return (mst_edges, sh_paths)


# def __pq_prim(adj_mat, u, edge_pq, visited, distances, incid_vec):
#    try:
#        __prim_scan(edge_pq, adj_mat, u, visited)
#        while True:
#            e = heapq.heappop(edge_pq)
#            u = e.head()
#            v = e.tail()
#
#            if visited[v] == 0:
#                incid_vec[v] = u
#                distances[v] = distances[u] + e.weight()
#                __prim_scan(edge_pq, adj_mat, v, visited)
#    except IndexError:
#        pass
#
#    return (incid_vec, distances)


def adj_mat_to_dijkstra_spath_mat(g):
    n = g.v()
    sh_path_mat = np.zeros((n, n))
    for i in range(n):
        mst_result = min_spanning_tree(g, i)
        for j in range(n):
            sh_path_mat[i, j] = mst_result[1][j]

    return sh_path_mat


def async_adj_mat_to_dijkstra_spath_mat(g):
    n = g.v()
    sh_path_mat = np.zeros((n, n))

    #    if os.path.exists('tmp.dat'):
    #        os.remove('tmp.dat')
    #    f = open('tmp.dat', 'wb')
    #    pickle.dump(g, f)
    #    f.close()

    processor_pool = mproc.Pool(mproc.cpu_count())
    async_results = [
        processor_pool.apply_async(min_spanning_tree, (g, i)) for i in range(n)
    ]
    processor_pool.close()
    processor_pool.join()
    for i in range(n):
        mst_result = async_results[i].get()
        for j in range(n):
            sh_path_mat[i, j] = mst_result[1][j]

    return sh_path_mat


def variance_weighting(sh_paths, pos_lab_indices, neg_lab_indices):
    pos_indices = []
    neg_indices = []
    pos_weights = []
    neg_weights = []

    for i in range(sh_paths.shape[0]):
        dist_to_pos = float(
            np.sum(sh_paths[i, pos_lab_indices]) / pos_lab_indices.shape[0]
        )
        dist_to_neg = float(
            np.sum(sh_paths[i, neg_lab_indices]) / neg_lab_indices.shape[0]
        )

        if dist_to_pos < dist_to_neg:
            pos_indices.append(i)
            pos_weights.append(dist_to_pos)
        else:
            neg_indices.append(i)
            neg_weights.append(dist_to_neg)

    pos_indices = np.array(pos_indices)
    neg_indices = np.array(neg_indices)
    pos_weights = np.array(pos_weights)
    neg_weights = np.array(neg_weights)
    pos_weights *= float(1.0 / np.sum(pos_weights))
    neg_weights *= float(1.0 / np.sum(neg_weights))

    return {"1": (pos_indices, pos_weights), "-1": (neg_indices, neg_weights)}


def largest_mst_weigth(data, incid_vec):
    sigma = 0.0
    for i in range(1, incid_vec.shape[0]):
        d = np.linalg.norm(data[i, :] - data[incid_vec[i], :])
        if d > sigma:
            sigma = d

    return sigma


def mean_pairwise_distance(X):
    euclidean_distances = dist.squareform(dist.pdist(X, "euclidean"))
    return np.mean(euclidean_distances)


def median_pairwise_distances(X):
    euclidean_distances = dist.squareform(dist.pdist(X, "euclidean"))
    return np.median(euclidean_distances)


def std_of_pairwise_distances(X):
    euclidean_distances = dist.squareform(dist.pdist(X, "euclidean"))
    return np.std(euclidean_distances)


def top_half_median(X, threshold):
    n = X.shape[0]
    th_pw_dists = []
    for i in range(n):
        for j in range(i):
            i2j_dist = np.linalg.norm(X[i, :] - X[j, :])
            if i2j_dist < threshold:
                th_pw_dists.append(i2j_dist)

    return np.median(th_pw_dists)


def mean_data_sq_norm(X):
    n = X.shape[0]
    tot_norm = 0.0
    for i in range(n):
        tot_norm = tot_norm + np.dot(X[i, :].reshape(1, -1), X[i, :].reshape(-1, 1))

    return float(tot_norm / n)


def __choose_marker(idx, idx_coll, label):
    for i in range(idx_coll.shape[0]):
        if idx == idx_coll[i]:
            return "s"

    if label == 1:
        return "o"

    return "x"


def planar_class_embedding(X, y, labelled_idxs, path):
    n = X.shape[0]
    d = X.shape[1]
    for i in range(n):
        x_coord = X[i, 0]
        y_coord = X[i, 1] if d > 1 else x_coord
        plt.plot(
            x_coord,
            y_coord,
            marker=__choose_marker(i, labelled_idxs, y[i, 0]),
            color="r" if y[i, 0] == 1 else "b",
            markersize=3,
        )
    plt.savefig(path, dpi=200)
    plt.close()
