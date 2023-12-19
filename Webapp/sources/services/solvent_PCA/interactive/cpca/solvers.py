"""
Created on Jan 23, 2014

@author: doglic
"""
import multiprocessing as mproc
import os
import pickle
import warnings

import numpy as np
from scipy.linalg import get_blas_funcs

from . import utils

# from scipy import linalg as spla
"""
    ############################### Asynchronous execution ###############################
"""


def async_rom_eig_vals(eig_sec_solver, args):
    try:
        return eig_sec_solver.root(*args)
    except:
        warnings.warn(
            "error while computing the eigenvalue using secular eigen-solver..."
        )
        return None


def async_rom_eig_vec(i, delta_vec, v, diag):
    try:
        t = diag[i] - diag
        t[i] = -1.0
        return np.sign(v[i]) * np.sqrt(np.abs(np.prod(delta_vec / t))) / delta_vec
    except:
        warnings.warn(
            "error while computing the eigenvector using secular-eigen solver..."
        )
        return None


"""
    ###############################                        ###############################
"""


class secular_solver(object):
    def __init__(self, precision, max_iters):
        self.precision = precision
        self.max_iters = max_iters

    def secular_eval(self, x, v, sigmas, r=1.0):
        pass

    def dx_secular_eval(self, x, v, sigmas):
        pass

    def root(self, i, v, sigmas, r=1.0):
        pass


class gander_secular_solver(secular_solver):
    def __init__(self, precision, max_iters):
        super(gander_secular_solver, self).__init__(precision, max_iters)
        self._default_eta = self.precision + 1e-8

    def dx_secular_eval(self, x, v, sigmas):
        return float(np.sum(2 * (np.power(v, 2) / np.power(sigmas - x, 3))))

    def secular_eval(self, x, v, sigmas, r):
        return np.sum(np.power(v / (sigmas - x), 2)) - r * r

    def _outer_shift(self, v, sigmas, r):
        eta = self._default_eta
        for i in range(v.shape[0]):
            if np.abs(v[-i - 1]) > 0.0:
                eta = np.abs(float(v[-i - 1] / r))
                break
        return eta

    def __is_outer_root(self, i, n):
        return i == -1 or i == n - 1

    def root(self, i, d, sigmas, r):
        if self.__is_outer_root(i, d.shape[0]):
            asc_sorted_sigmas = np.argsort(sigmas)
            root = self._outer_shift(d[asc_sorted_sigmas], sigmas[asc_sorted_sigmas], r)
            squared_radius = r * r
            shifted_sigmas = sigmas - sigmas[asc_sorted_sigmas[-1]]
            eta = 1.0
            k = 0
            while k < self.max_iters:
                u = self.secular_eval(root, d, shifted_sigmas, r) + squared_radius
                v = self.dx_secular_eval(root, d, shifted_sigmas)
                eta = 2 * float(u / v) * (1.0 - float(np.sqrt(u) / r))

                if eta < self.precision:
                    break
                elif np.abs(root) < self.precision:
                    warnings.warn("gander secular solver converged to the end point!")
                    break
                elif np.abs(v) < self.precision:
                    warnings.warn(
                        "gander secular solver has derivative 0 at the candidate point!"
                    )
                    break

                root += eta
                k += 1

            if k == self.max_iters:
                at_root = self.secular_eval(root, d, shifted_sigmas, r)
                warnings.warn(
                    "max number of allowed iterations reached with the gander secular solver (value at root = "
                    + str(at_root)
                    + ")"
                )

            return (root + sigmas[asc_sorted_sigmas[-1]], shifted_sigmas - root)

        raise ValueError(
            "Operation not supported! It is possible to compute only the largest root!"
        )


class eigen_secular_solver(secular_solver):
    def __init__(self, precision, max_iters):
        super(eigen_secular_solver, self).__init__(precision, max_iters)
        self._default_outer_shift = self.precision + 1e-8

    def _phi_eval(self, x, v, sigmas):
        if v.shape[0] == 0:
            return 0.0
        return np.sum((v * v) / (sigmas - x))

    def _dx_phi_eval(self, x, v, sigmas):
        if v.shape[0] == 0:
            return 0.0
        return np.sum(np.power(v / (sigmas - x), 2))

    def secular_eval(self, x, v, sigmas, r=1.0):
        return float(1.0 + self._phi_eval(x, v, sigmas))

    def dx_secular_eval(self, x, v, sigmas):
        return self._dx_phi_eval(x, v, sigmas)

    def __is_outer_root(self, i, n):
        return i == -1 or i == (n - 1)

    def _outer_shift(self, v):
        for i in range(v.shape[0]):
            if np.abs(v[-1 - i]) > 0.0:
                return v[-1 - i] * v[-1 - i]
        return self.precision

    def _eta(self, a, b, c):
        if a > 0:
            return float((2.0 * b) / (a + np.sqrt(a * a - 4 * b * c)))
        else:
            return float((a - np.sqrt(a * a - 4 * b * c)) / (2.0 * c))

    def _inner_init_root(self, v, sigmas, split_idx, axis_shift):
        midpoint = float((sigmas[split_idx - 1] + sigmas[split_idx]) / 2.0)
        c = self.secular_eval(
            midpoint, v[: (split_idx - 1)], sigmas[: (split_idx - 1)]
        ) + self._phi_eval(midpoint, v[(split_idx + 1) :], sigmas[(split_idx + 1) :])
        delta0 = sigmas[split_idx] - sigmas[split_idx - 1]
        if np.abs(axis_shift - sigmas[split_idx]) < self.precision:
            a = (
                -c * delta0
                + v[split_idx - 1] * v[split_idx - 1]
                + v[split_idx] * v[split_idx]
            )
            b = -v[split_idx] * v[split_idx] * delta0
        else:
            a = (
                c * delta0
                + v[split_idx - 1] * v[split_idx - 1]
                + v[split_idx] * v[split_idx]
            )
            b = v[split_idx - 1] * v[split_idx - 1] * delta0

        return self._eta(a, b, c)

    def _inner_root(self, v, sigmas, split_idx):
        eval_precision = sigmas.shape[0] * self.precision
        lep = sigmas[split_idx - 1]
        rep = sigmas[split_idx]
        if np.abs(lep - rep) < self.precision:
            raise ValueError("diagonal matrix must be deflated!")

        root = float((lep + rep) / 2.0)
        eta = 1.0
        axis_shift = 0.0
        at_midpoint = self.secular_eval(root, v, sigmas)
        if at_midpoint > eval_precision:
            axis_shift = lep
        elif at_midpoint < -eval_precision:
            axis_shift = rep
        else:
            return root

        root = self._inner_init_root(v, sigmas, split_idx, axis_shift)
        rep -= axis_shift
        lep -= axis_shift
        shifted_sigmas = sigmas - axis_shift
        lep_delta = lep - root
        rep_delta = rep - root
        acceptance_threshold = eval_precision
        sec_at_root = 1.0
        itr = 0
        while itr < self.max_iters:
            lphi = self._phi_eval(root, v[:split_idx], shifted_sigmas[:split_idx])
            dx_lphi = self._dx_phi_eval(root, v[:split_idx], shifted_sigmas[:split_idx])
            rphi = self._phi_eval(root, v[split_idx:], shifted_sigmas[split_idx:])
            dx_rphi = self._dx_phi_eval(root, v[split_idx:], shifted_sigmas[split_idx:])
            sec_at_root = 1.0 + rphi + lphi
            if np.abs(sec_at_root) < acceptance_threshold:
                break

            sec_prec_scaler = 1.0 + np.abs(rphi) + np.abs(lphi)
            delta_prod = lep_delta * rep_delta
            a = (lep_delta + rep_delta) * sec_at_root - (dx_rphi + dx_lphi) * delta_prod
            b = sec_at_root * delta_prod
            c = sec_at_root - lep_delta * dx_lphi - rep_delta * dx_rphi
            eta = self._eta(a, b, c)
            root += eta
            lep_delta = float(lep - root)
            rep_delta = float(rep - root)
            acceptance_threshold = max(eval_precision * sec_prec_scaler, self.precision)

            # !!! bounce back into the interval !!!
            if lep_delta * rep_delta > 0.0:
                warnings.warn("bouncing the root back into the interval...")
                if lep_delta > 0:
                    root += -eta + 0.5 * np.abs(rep_delta + eta)
                else:
                    root += -eta + 0.5 * np.abs(lep_delta + eta)
            elif np.abs(lep_delta) < self.precision:
                #                warnings.warn(str(axis_shift) + ' - closing to the left end point of the interval (' + str(sigmas[split_idx - 1]) + ', ' +  str(sigmas[split_idx]) + ')! Delta = ' + str(lep_delta))
                break
            elif np.abs(rep_delta) < self.precision:
                #                warnings.warn(str(axis_shift) + ' - closing to the right end point of the interval (' + str(sigmas[split_idx - 1]) + ', ' +  str(sigmas[split_idx]) + ')! Delta = ' + str(rep_delta))
                break

            itr += 1

        if itr == self.max_iters:
            warnings.warn(
                "max number of allowed iterations reached with the eigen secular solver (value at _inner_ root = "
                + str(sec_at_root)
                + ")"
            )

        return (root + axis_shift, shifted_sigmas - root)

    def _outer_root(self, v, sigmas):
        shifted_sigmas = sigmas - sigmas[-1]
        root = self._outer_shift(v)
        sec_at_root = 1.0
        eta = 1.0
        itr = 0
        while itr < self.max_iters:
            phi = self._phi_eval(root, v, shifted_sigmas)
            dx_phi = self._dx_phi_eval(root, v, shifted_sigmas)
            sec_at_root = phi + 1.0
            eta = (phi * sec_at_root) / dx_phi
            if eta < self.precision:
                break
            root += eta
            itr += 1

        if itr == self.max_iters:
            warnings.warn(
                "max number of allowed iterations reached with the eigen secular solver (value at _outer_ root = "
                + str(sec_at_root)
                + ")"
            )

        return (root + sigmas[-1], shifted_sigmas - root)

    def root(self, i, v, sigmas, r=1.0):
        n = sigmas.shape[0]
        if self.__is_outer_root(i, n):
            return self._outer_root(v, sigmas)
        return self._inner_root(v, sigmas, i + 1)


# class graggs_eigen_secular_solver(eigen_secular_solver):
#
#    def __init__(self, precision, max_iters):
#        super(graggs_eigen_secular_solver, self).__init__(precision, max_iters)
#
#    def dx_2_secular_eval(self, x, v, sigmas):
#        return 2.0 * np.sum(v * v / np.power(sigmas - x, 3))
#
#    def _inner_root(self, v, sigmas, split_idx):
#        eval_precision = sigmas.shape[0] * self.precision
#        lep = sigmas[split_idx - 1]
#        rep = sigmas[split_idx]
#        if np.abs(lep - rep) < self.precision:
#            raise ValueError, 'diagonal matrix must be deflated!'
#
#        root = float((lep + rep) / 2.0)
#        eta = 1.0
#        axis_shift = 0.0
#        at_axis_shift_decider = self.secular_eval(root, v, sigmas)
#        c0 = at_axis_shift_decider - self._phi_eval(root, v[(split_idx-1):(split_idx+1)], sigmas[(split_idx-1):(split_idx+1)])
#        delta0 = rep - lep
#        if at_axis_shift_decider > eval_precision:
#            axis_shift = lep
#            a0 = c0 * delta0 + v[split_idx - 1] * v[split_idx - 1] + v[split_idx] * v[split_idx]
#            b0 = v[split_idx - 1] * v[split_idx - 1] * delta0
#        elif at_axis_shift_decider < -eval_precision:
#            axis_shift = rep
#            a0 = -c0 * delta0 + v[split_idx - 1] * v[split_idx - 1] + v[split_idx] * v[split_idx]
#            b0 = -v[split_idx] * v[split_idx] * delta0
#        else:
#            return root
#
#        if a0 > 0:
#            root = float((2.0 * b0) / (a0 + np.sqrt(a0 * a0 - 4 * b0 * c0)))
#        else:
#            root = float((a0 - np.sqrt(a0 * a0 - 4 * b0 * c0)) / (2.0 * c0))
#
#        rep -= axis_shift
#        lep -= axis_shift
#        shifted_sigmas = sigmas - axis_shift
#        lep_delta = lep - root
#        rep_delta = rep - root
#        acceptance_threshold = eval_precision
#        sec_at_root = 1.0
#        itr = 0
#        while np.abs(sec_at_root) > acceptance_threshold and itr < self.max_iters:
#            dx_sec_at_root = self.dx_secular_eval(root, v, shifted_sigmas)
#            dx_2_at_root = self.dx_2_secular_eval(root, v, shifted_sigmas)
#            lphi = self._phi_eval(root, v[:split_idx], shifted_sigmas[:split_idx])
#            rphi = self._phi_eval(root, v[split_idx:], shifted_sigmas[split_idx:])
#            s = ((np.power(lep_delta, 3) * rep_delta) / (lep_delta - rep_delta)) * (dx_sec_at_root / rep_delta - 0.5 * dx_2_at_root)
#            S = ((lep_delta * np.power(rep_delta, 3)) / (rep_delta - lep_delta)) * (dx_sec_at_root / lep_delta - 0.5 * dx_2_at_root)
#            sec_at_root = 1.0 + rphi + lphi
#            sec_prec_scaler = 1.0 + np.abs(rphi) + np.abs(lphi)
#            delta_prod = lep_delta * rep_delta
#
#            c = sec_at_root - (lep_delta + rep_delta) * dx_sec_at_root + 0.5 * delta_prod * dx_2_at_root
#            a = c * (lep_delta + rep_delta) + s + S
#            b = c * delta_prod + s * rep_delta + S * lep_delta
#
#            if a > 0:
#                eta = float((2 * b) / (a + np.sqrt(a * a - 4 * b * c)))
#            else:
#                eta = float((a - np.sqrt(a * a - 4 * b * c)) / (2 * c))
#
#            root += eta
#            lep_delta = float(lep - root)
#            rep_delta = float(rep - root)
#            acceptance_threshold = max(eval_precision * sec_prec_scaler, self.precision)
#
#            # !!! bounce back into the interval !!!
#            if lep_delta * rep_delta > 0.0:
#                warnings.warn('bouncing the root back into the interval...')
#                if lep_delta > 0:
#                    root += -eta + 0.5 * rep_delta
#                else:
#                    root += -eta + 0.5 * lep_delta
#            elif np.abs(lep_delta) < self.precision:
#                warnings.warn(str(axis_shift) + ' - closing to the left end point of the interval (' + str(sigmas[split_idx - 1]) + ', ' +  str(sigmas[split_idx]) + ')! Delta = ' + str(lep_delta))
#                root += -eta + float(0.5 * np.abs(rep_delta + eta)) + axis_shift - sigmas[split_idx]
#                lep += axis_shift - sigmas[split_idx]
#                rep += axis_shift - sigmas[split_idx]
#                lep_delta = float(lep - root)
#                rep_delta = float(rep - root)
#                shifted_sigmas += axis_shift - sigmas[split_idx]
#                axis_shift = sigmas[split_idx]
#            elif np.abs(rep_delta) < self.precision:
#                warnings.warn(str(axis_shift) + ' - closing to the right end point of the interval (' + str(sigmas[split_idx - 1]) + ', ' +  str(sigmas[split_idx]) + ')! Delta = ' + str(rep_delta))
#                root += -eta - float(0.5 * np.abs(lep_delta + eta)) + axis_shift - sigmas[split_idx - 1]
#                lep += axis_shift - sigmas[split_idx - 1]
#                rep += axis_shift - sigmas[split_idx - 1]
#                lep_delta = float(lep - root)
#                rep_delta = float(rep - root)
#                shifted_sigmas += axis_shift - sigmas[split_idx - 1]
#                axis_shift = sigmas[split_idx - 1]
#
#            itr += 1
#
#        if itr == self.max_iters:
#            warnings.warn('max number of allowed iterations reached with the eigen secular solver (value at _inner_ root = ' + str(sec_at_root) + ')')
#
#        return root + axis_shift


class givens_transformer(object):
    def __init__(self, precision=2.56e-16):
        self.__precision = precision

    def compute_coeffs(self, a, b):
        c = 0.0
        s = 0.0
        if np.abs(b) < self.__precision:
            c = float(np.sign(a))
        elif np.abs(b) > np.abs(a):
            tau = -float(a / b)
            s = -float(1.0 / np.sqrt(1 + tau * tau)) * float(np.sign(b))
            c = s * tau
        else:
            tau = -float(b / a)
            c = float(1.0 / np.sqrt(1 + tau * tau)) * float(np.sign(a))
            s = c * tau

        return (c, s)

    def transform(self, A, i, j, givens_coeff):
        rotation = np.array(
            [[givens_coeff[0], -givens_coeff[1]], [givens_coeff[1], givens_coeff[0]]]
        )
        A[[i, j], :] = rotation.dot(A[[i, j], :])
        return A


class appended_column_qr_solver(object):
    def __init__(self, precision=2.56e-15):
        self.__givens_transf = givens_transformer()
        self.__precision = precision

    def decompose(self, A, Q, R, k):
        n = Q.shape[0]
        R_ = np.copy(R)
        Q_ = np.copy(Q).T
        v = Q_.dot(A[:, k].reshape(-1, 1))
        for i in reversed(range(k + 1, n)):
            givens_coeff = self.__givens_transf.compute_coeffs(v[i - 1, 0], v[i, 0])
            v = self.__givens_transf.transform(v, i - 1, i, givens_coeff)
            Q_ = self.__givens_transf.transform(Q_, i - 1, i, givens_coeff)
        R_[: (k + 1), k] = v[: (k + 1), 0]
        return (Q_.T, R_)

    def qr_transform(self, W, Q, R, alpha):
        n = Q.shape[0]
        m = 1 + R.shape[1] if R is not None else 1
        R_ = np.zeros((n, m))
        R_[:, : (m - 1)] = R
        Q = Q.T
        v = Q.dot(alpha.reshape(-1, 1))
        for i in reversed(range(m, n)):
            givens_coeff = self.__givens_transf.compute_coeffs(v[i - 1, 0], v[i, 0])
            v = self.__givens_transf.transform(v, i - 1, i, givens_coeff)
            Q = self.__givens_transf.transform(Q, i - 1, i, givens_coeff)
            W = self.__givens_transf.transform(W, i - 1, i, givens_coeff)
        R_[:m, m - 1] = v[:m, 0]
        return (Q.T, R_, W)


class rom_des_eigen_solver(object):
    def __init__(self, precision, max_iters):
        self._eig_solver = eigen_secular_solver(precision, max_iters)
        #        self._eig_solver = graggs_eigen_secular_solver(precision, max_iters)
        self.__givens_transformer = givens_transformer()
        self.precision = precision
        self.gu_eisen_threshold = 1e2

    def __inv_arg_permutation(self, arg_perm, n):
        inv_arg_perm = list(range(n))
        for i in range(n):
            inv_arg_perm[arg_perm[i]] = i
        return inv_arg_perm

    def _eig_vals(self, v, diag):
        n = v.shape[0]
        eig_vals = np.zeros(n)
        delta_vecs = np.zeros((n, n))
        for i in range(n):
            res = self._eig_solver.root(i, v, diag)
            eig_vals[i] = res[0]
            delta_vecs[:, i] = res[1]
        return (eig_vals, delta_vecs)

    def _eig_vecs(self, delta_vecs, v, diag):
        n = v.shape[0]
        eig_vecs = np.ones((n, n))
        for i in range(n):
            t = diag[i] - diag
            t[i] = -1.0
            eig_vecs[i, :] = (
                np.sign(v[i])
                * np.sqrt(np.abs(np.prod(delta_vecs[i, :] / t)))
                / delta_vecs[i, :]
            )
        #            eig_vecs[i, :] = permuted_v[i] / s

        for i in range(n):
            eig_vecs[:, i] = eig_vecs[:, i] / np.linalg.norm(eig_vecs[:, i])

        return eig_vecs

    def __extract_deflation_blocks(self, diag, deflate_precision):
        givens_deflates = []
        block_repres = None
        block_deflator = []
        for i in range(diag.shape[0]):
            if (
                block_repres is None
                or np.abs(diag[i] - block_repres) > deflate_precision
            ):
                block_repres = diag[i]
                if len(block_deflator) > 1:
                    givens_deflates.append(block_deflator)
                block_deflator = [i]
            else:
                block_deflator.append(i)
        if len(block_deflator) > 1:
            givens_deflates.append(block_deflator)

        return givens_deflates

    def __inv_deflate_eig_vecs(self, eig_vecs, givens_deflates, deflate_precision):
        for gd in reversed(givens_deflates):
            givens_coeff = (gd[0][0], -gd[0][1])
            eig_vecs = self.__givens_transformer.transform(
                eig_vecs, gd[1], gd[2], givens_coeff
            )
        return eig_vecs

    def __deflate(self, v, diag, deflate_precision):
        n = diag.shape[0]
        v_ = np.copy(v)
        deflation_blocks = self.__extract_deflation_blocks(diag, deflate_precision)
        nulled_poss = []
        givens_deflates = []

        for db in deflation_blocks:
            m = len(db)
            for i in reversed(range(1, m)):
                if np.abs(v_[db[i]]) < deflate_precision:
                    nulled_poss.append(db[i])
                    v_[db[i]] = 0.0
                else:
                    givens_coeff = self.__givens_transformer.compute_coeffs(
                        v_[db[i - 1]], v_[db[i]]
                    )
                    v_ = self.__givens_transformer.transform(
                        v_.reshape(-1, 1), db[i - 1], db[i], givens_coeff
                    )
                    nulled_poss.append(db[i])
                    givens_deflates.append((givens_coeff, db[i - 1], db[i]))
                    v_[
                        db[i]
                    ] = 0.0  # it is below the threshold due to the Givens rotation

        nnp = np.sort(list(set(range(n)) - set(nulled_poss)))

        non_nulled_pos = []
        for i in nnp:
            if np.abs(v_[i]) < deflate_precision:
                v_[i] = 0.0
                nulled_poss.append(i)
            else:
                non_nulled_pos.append(i)

        return (v_, non_nulled_pos, nulled_poss, givens_deflates)

    def decompose(self, diag, v, ro=1.0, sorted_eig_sys=False):
        n = diag.shape[0]
        sign_ro = float(np.sign(ro))
        v_scaler = np.sqrt(np.abs(ro))
        asc_diag_args = np.argsort(sign_ro * diag)
        permuted_diag = sign_ro * diag[asc_diag_args]
        permuted_v = v_scaler * v[asc_diag_args]

        cond_num = float(
            np.linalg.norm(permuted_diag)
            + permuted_v.reshape(1, -1).dot(permuted_v.reshape(-1, 1))
        )
        deflate_precision = self.precision * cond_num
        deflation_res = self.__deflate(permuted_v, permuted_diag, deflate_precision)
        nontrivial_v = (deflation_res[0][deflation_res[1]]).reshape(-1)
        nontrivial_diag = (permuted_diag[deflation_res[1]]).reshape(-1)

        eig_val_res = self._eig_vals(nontrivial_v, nontrivial_diag)
        deflated_eig_vecs = self._eig_vecs(
            eig_val_res[1], nontrivial_v, nontrivial_diag
        )
        eig_vals = np.zeros(n)
        eig_vals[deflation_res[1]] = eig_val_res[0]
        eig_vals[deflation_res[2]] = permuted_diag[deflation_res[2]]
        eig_vecs = np.identity(n)
        for i, c in enumerate(deflation_res[1]):
            eig_vecs[deflation_res[1], c] = deflated_eig_vecs[:, i]
        eig_vecs = self.__inv_deflate_eig_vecs(
            eig_vecs, deflation_res[3], deflate_precision
        )

        if not sorted_eig_sys:
            inv_arg_perm = self.__inv_arg_permutation(asc_diag_args, n)
            eig_vals = eig_vals[inv_arg_perm]
            eig_vecs = eig_vecs[:, inv_arg_perm][inv_arg_perm, :]
        eig_vals = sign_ro * eig_vals

        return (eig_vals, eig_vecs)


class async_rom_des_eigen_solver(rom_des_eigen_solver):
    def __init__(self, precision, max_iters, num_threads):
        super(async_rom_des_eigen_solver, self).__init__(precision, max_iters)
        self._num_threads = num_threads

    def _eig_vals(self, v, diag):
        n = v.shape[0]
        eig_vals = np.zeros(n)
        delta_vecs = np.zeros((n, n))

        processor_pool = mproc.Pool(self._num_threads)
        async_results = [
            processor_pool.apply_async(
                async_rom_eig_vals, (self._eig_solver, (i, v, diag))
            )
            for i in range(n)
        ]
        processor_pool.close()
        processor_pool.join()
        for i in range(n):
            res = async_results[i].get()
            eig_vals[i] = res[0]
            delta_vecs[:, i] = res[1]

        return (eig_vals, delta_vecs)

    def _eig_vecs(self, delta_vecs, v, diag):
        n = v.shape[0]
        eig_vecs = np.ones((n, n))

        processor_pool = mproc.Pool(self._num_threads)
        async_results = [
            processor_pool.apply_async(
                async_rom_eig_vec, (i, delta_vecs[i, :], v, diag)
            )
            for i in range(n)
        ]
        processor_pool.close()
        processor_pool.join()
        for i in range(n):
            eig_vecs[i, :] = async_results[i].get()

        for i in range(n):
            eig_vecs[:, i] = eig_vecs[:, i] / np.linalg.norm(eig_vecs[:, i])

        return eig_vecs


class quad_solver(object):
    def __init__(self, precision, max_iters):
        self.gander_solver = gander_secular_solver(precision, max_iters)
        self.precision = precision
        self.__degenerate_sol_precision = 1e-5

    def __inv_arg_permutation(self, arg_perm, n):
        inv_arg_perm = range(n)
        for i in range(n):
            inv_arg_perm[arg_perm[i]] = i
        return inv_arg_perm

    def maximizer(self, eig_sys, b, r):
        top_eig_val_idx = np.argmax(eig_sys[0])
        if b is not None:
            d = eig_sys[1].T.dot(b.reshape(-1, 1)).reshape(-1)
            gander_res = self.gander_solver.root(-1, d, eig_sys[0], r)
            if np.abs(gander_res[1][top_eig_val_idx]) > self.precision:
                return (
                    eig_sys[1].dot((d / gander_res[1]).reshape(-1, 1)).reshape(-1),
                    d,
                )

            diag = np.copy(gander_res[1])
            diag[:top_eig_val_idx] = 1.0 / diag[:top_eig_val_idx]
            diag[top_eig_val_idx] = 0.0
            diag[(top_eig_val_idx + 1) :] = 1.0 / diag[(top_eig_val_idx + 1) :]

            pinv_solution = (
                eig_sys[1].dot((d.reshape(-1) * diag).reshape(-1, 1)).reshape(-1)
            )
            pinv_sol_norm = np.linalg.norm(pinv_solution)
            if np.abs(pinv_sol_norm - r) < self.__degenerate_sol_precision:
                warnings.warn(
                    "unique solution (using pseudo-inverse) computed successfully!"
                )
                return (pinv_solution, d)
            elif pinv_sol_norm < r:
                warnings.warn("a non-unique solution computed with the help of pinv!")
                norm_delta = r - pinv_sol_norm
                pinv_solution[top_eig_val_idx] += norm_delta
                return (pinv_solution, d)

            raise ValueError(
                "gander lambda does not exist - the optimization problem is infeasible!"
            )

        return (r * eig_sys[1][:, top_eig_val_idx], None)


class greedy_dir_solver(object):
    def __init__(self, precision, max_iters, parallelize):
        self.quad_solver = quad_solver(precision, max_iters)
        self.qr_solver = appended_column_qr_solver(precision)
        self.precision = precision
        if parallelize:
            self.eig_solver = async_rom_des_eigen_solver(
                precision, max_iters, mproc.cpu_count()
            )
        else:
            self.eig_solver = rom_des_eigen_solver(precision, max_iters)

    def _next_direction(self, quad_eig_sys, b, r, alpha=None, mu=None):
        if alpha is None:
            quad_res = self.quad_solver.maximizer(quad_eig_sys, b, r)
            return (quad_res[0], quad_res[1], quad_eig_sys)

        rnk_one_upd_eig_sys = self.eig_solver.decompose(quad_eig_sys[0], alpha, mu)
        quad_res = self.quad_solver.maximizer(rnk_one_upd_eig_sys, b, r)
        return (quad_res[0], quad_res[1], rnk_one_upd_eig_sys)

    def from_sphere_to_ellipsoid(self, alpha, kernel_sys):
        if kernel_sys is not None:
            return kernel_sys[2].dot(alpha)
        return alpha

    def __nonzero_vec(self, v):
        for i in range(v.shape[0]):
            if np.abs(v[i]) > self.precision:
                return True
        return False

    def directions(self, kernel_sys, quad_eig_sys, b, r, orth_nu):
        n = quad_eig_sys[0].shape[0]
        dim = b.shape[1]
        directions = np.zeros((n, dim))
        d = np.copy(b)
        if isinstance(quad_eig_sys[1], list):
            dir_eig_sys = [(quad_eig_sys[0], quad_eig_sys[1][-1])]
            for i in range(len(quad_eig_sys[1]) - 1):
                d = quad_eig_sys[1][i].T.dot(d)
        else:
            dir_eig_sys = [quad_eig_sys]

        alpha = None
        for i in range(dim):
            dir_res = self._next_direction(
                dir_eig_sys[i - 1],
                d[:, i] if self.__nonzero_vec(d[:, i]) else None,
                r,
                alpha,
                orth_nu,
            )
            directions[:, i] = dir_res[0]
            for j in reversed(range(1, i + 1)):
                directions[:, i] = dir_eig_sys[j][1].dot(directions[:, i])

            if isinstance(quad_eig_sys[1], list):
                for j in reversed(range(len(quad_eig_sys[1]) - 1)):
                    directions[:, i] = quad_eig_sys[1][j].dot(directions[:, i])

            directions[:, i] = self.from_sphere_to_ellipsoid(
                directions[:, i], kernel_sys
            )

            d = dir_res[2][1].T.dot(d)
            alpha = np.copy(dir_res[0])
            alpha = dir_res[2][1].T.dot(alpha)
            dir_eig_sys.append(dir_res[2])

        return directions

    def _seq_eig_sys_next_direction(self, seq_eig_sys, b, r):
        k = len(seq_eig_sys[1])
        d = np.copy(b)
        for i in range(k - 1):
            d = seq_eig_sys[1][i].T.dot(d.reshape(-1, 1))

        sol = self.quad_solver.maximizer((seq_eig_sys[0], seq_eig_sys[1][-1]), d, r)[0]
        for i in reversed(range(k - 1)):
            sol = seq_eig_sys[1][i].dot(sol.reshape(-1, 1))

        return sol.reshape(-1)

    def __horth_qr_transform(self, W, b, m, qr_res):
        ort_map_W = qr_res[0].T.dot(W).dot(qr_res[0])
        C = ort_map_W[m:, m:]
        f = qr_res[0].T.dot(b)[m:, :].reshape(-1, 1)
        return (C, f)

    def _next_horth_direction(self, eig_sys, b, r, qr_sys, alpha):
        if None is qr_sys and None is alpha:
            return (
                self._next_direction(eig_sys, b if self.__nonzero_vec(b) else None, r)[
                    0
                ],
                None,
                eig_sys,
            )

        if qr_sys is None:
            #            Q, R = spla.qr(alpha.reshape(-1, 1))
            #            qr_sys = (Q, R)
            #            U = Q.T.dot(eig_sys[1])
            res = self.qr_solver.qr_transform(
                eig_sys[1], np.identity(alpha.shape[0]), None, alpha
            )
            qr_sys = (res[0], res[1])
            U = res[2]
        else:
            res = self.qr_solver.qr_transform(eig_sys[1], qr_sys[0], qr_sys[1], alpha)
            qr_sys = (res[0], res[1])
            U = res[2]

        m = qr_sys[1].shape[1]
        eig_sys = (eig_sys[0], U)
        d = qr_sys[0].T.dot(b.reshape(-1, 1)).reshape(-1)[m:]
        u22 = eig_sys[1][m:, :][:, m:]
        u21 = eig_sys[1][m:, :m] if m > 1 else eig_sys[1][m:, :m].reshape(-1, 1)
        diag = eig_sys[0][m:]
        eig_vecs = [u22]
        for i in range(m):
            v = np.copy(u21[:, i]).reshape(-1, 1)
            for ev in eig_vecs:
                v = ev.T.dot(v)
            res_eig_sys = self.eig_solver.decompose(diag, v.reshape(-1), eig_sys[0][i])
            diag = res_eig_sys[0]
            eig_vecs.append(res_eig_sys[1])

        sub_sol = self._seq_eig_sys_next_direction((diag, eig_vecs), d, r)
        #        evl, evc = np.linalg.eigh(U.dot(np.diag(eig_sys[0])).dot(U.T)[m:, :][:, m:])
        #        sub_sol = self.quad_solver.maximizer((evl, evc), d, r)[0]

        sol = np.zeros(alpha.shape[0])
        sol[m:] = sub_sol
        sol = qr_sys[0].dot(sol.reshape(-1, 1)).reshape(-1)
        return (sol, qr_sys, eig_sys)

    def hard_orth_directions(self, kernel_sys, quad_eig_sys, b, r):
        n = quad_eig_sys[0].shape[0]
        dim = b.shape[1]
        directions = np.zeros((n, dim))
        d = np.copy(b)
        qr_sys = None
        alpha = None
        eig_sys = (np.copy(quad_eig_sys[0]), np.copy(quad_eig_sys[1]))
        for i in range(dim):
            dir_res = self._next_horth_direction(eig_sys, d[:, i], r, qr_sys, alpha)
            alpha = np.copy(dir_res[0]).reshape(-1)
            qr_sys = dir_res[1]
            eig_sys = dir_res[2]
            directions[:, i] = self.from_sphere_to_ellipsoid(dir_res[0], kernel_sys)

        return directions


class embedder(object):
    def __init__(self, precision, max_iters, parallelize=False):
        self.direction_solver = greedy_dir_solver(precision, max_iters, parallelize)
        if parallelize:
            self.eig_solver = async_rom_des_eigen_solver(
                precision, max_iters, mproc.cpu_count()
            )
        else:
            self.eig_solver = rom_des_eigen_solver(precision, max_iters)

    def kernel_sys(self, K):
        eigvals, eigvecs = np.linalg.eigh(K)
        min_real = float(min(np.real(eigvals)))
        max_real = np.abs(max(np.real(eigvals)))
        mat_cond_num = max_real / np.abs(min_real)
        if mat_cond_num > 1e6:
            warnings.warn(
                "WARNING: The kernel matrix K is poorly conditioned --- cond_num = "
                + str(mat_cond_num)
                + ". Addding noise to the diagonal..."
            )
            eigvals = eigvals + 2.56e-8
        elif min_real <= 0:
            warnings.warn(
                "WARNING: Numerically unstable eigendecomposition of p.d. matrix K --- min_eigen_val = "
                + str(min_real)
                + ". Addding noise to diagonal..."
            )
            eigvals = eigvals + np.abs(min_real) + 2.56e-8

        sqrt_eig_vals = np.sqrt(eigvals)
        sqrt_inv_diag = 1.0 / sqrt_eig_vals
        sqrt_eig_vals = sqrt_eig_vals.reshape(1, -1)
        sqrt_inv_diag = sqrt_inv_diag.reshape(1, -1)

        sqrt_mat = (eigvecs * sqrt_eig_vals).dot(eigvecs.T)
        sqrt_inv_mat = (eigvecs * sqrt_inv_diag).dot(eigvecs.T)
        return (K, sqrt_mat, sqrt_inv_mat, eigvals, eigvecs)

    def sph_cl_var_term_eig_sys(self, kernel_sys):
        n = kernel_sys[0].shape[0]
        # H = np.identity(n) - float(1.0 / n)# * np.ones((n, n))
        # USE BLAS MATRIX MULTIPLICATION
        gemm = get_blas_funcs(["gemm"], [kernel_sys[0]])
        K2 = gemm[0](1, kernel_sys[0], kernel_sys[0])
        # K2 = kernel_sys[0].dot(kernel_sys[0])
        x = np.sum(kernel_sys[0], axis=0).reshape(-1, 1)
        X = ((1.0 / n) * x).dot(x.T)
        sph_var_term = K2 - X
        # sph_var_term = kernel_sys[1].dot(H).dot(kernel_sys[1])
        svt_eig_vals, svt_eig_vecs = np.linalg.eigh(sph_var_term)
        return (svt_eig_vals, svt_eig_vecs)

    def unl_sph_cl_var_term_eig_sys(self, kernel_sys, cp_indxs, mu):
        m = cp_indxs.shape[0]
        n = kernel_sys[0].shape[0]
        H = np.ones((n, n))
        H[:m, :m] = np.identity(m) + float((m - 2 * n) / (n * n)) * np.ones((m, m))
        H[:m, m:] = -float(1.0 / n) * np.ones((m, n - m))
        H[m:, :m] = -float(1.0 / n) * np.ones((n - m, m))
        H[m:, m:] = float(m / (n * n)) * np.ones((n - m, n - m))
        sph_var_term = kernel_sys[1].dot(H).dot(kernel_sys[1])
        tmp_mat = kernel_sys[2].dot(kernel_sys[0][:, cp_indxs])
        reg_block = mu * tmp_mat.dot(tmp_mat.T)
        svt_eig_vals, svt_eig_vecs = np.linalg.eigh(sph_var_term - reg_block)
        return (svt_eig_vals, svt_eig_vecs)

    def sph_cp_quad_term_eig_sys(self, kernel_sys, quad_eig_sys, new_cp_idx, mu):
        v = kernel_sys[2].dot(kernel_sys[0][:, new_cp_idx])
        eig_vec_seq = []
        if isinstance(quad_eig_sys[1], list):
            for tes in quad_eig_sys[1]:
                v = tes.T.dot(v)
            eig_vec_seq += quad_eig_sys[1]
        else:
            v = quad_eig_sys[1].T.dot(v)
            eig_vec_seq.append(quad_eig_sys[1])

        rou_eig_sys = self.eig_solver.decompose(quad_eig_sys[0], v.reshape(-1), -mu)

        eig_vec_seq.append(rou_eig_sys[1])
        #        U = quad_eig_sys[1].dot(rou_eig_sys[1]) # !!! avoid matrix multiplication
        return (rou_eig_sys[0], eig_vec_seq)

    def sph_cl_weighted_var_term_eig_sys(
        self, X, y, label_mask, kernel_sys, sh_paths_cache="sh_paths_cache.dat"
    ):
        if os.path.exists(sh_paths_cache):
            f = open(sh_paths_cache, "rb")
            sh_paths = pickle.load(f)
            f.close()
        else:
            edist_adj_mat = utils.pairwise_distance_mat(X)
            g = utils.sim_mat_to_graph(edist_adj_mat)
            mst_result = utils.min_spanning_tree(g, 0)
            cutoff = utils.largest_mst_weigth(X, mst_result[0]) + 1e-6
            g = utils.sim_mat_to_graph(edist_adj_mat, cutoff, False)
            sh_paths = utils.adj_mat_to_dijkstra_spath_mat(g)

            f = open(sh_paths_cache, "wb")
            pickle.dump(sh_paths, f)
            f.close()
        #        sh_paths = utils.async_adj_mat_to_dijkstra_spath_mat(g)

        pos_lab_indices = []
        neg_lab_indices = []
        for i in range(y.shape[0]):
            if y[i, 0] > 0:
                pos_lab_indices.append(label_mask[i])
            else:
                neg_lab_indices.append(label_mask[i])
        pos_lab_indices = np.array(pos_lab_indices)
        neg_lab_indices = np.array(neg_lab_indices)

        var_wdict = utils.variance_weighting(sh_paths, pos_lab_indices, neg_lab_indices)
        pos_dict_elem = var_wdict["1"]
        neg_dict_elem = var_wdict["-1"]

        pos_sub_mat_kernel = kernel_sys[0][pos_dict_elem[0], :]
        neg_sub_mat_kernel = kernel_sys[0][neg_dict_elem[0], :]
        pos_cl_weigths = np.diag(pos_dict_elem[1])
        neg_cl_weigths = np.diag(neg_dict_elem[1])

        n1 = pos_dict_elem[0].shape[0]
        n2 = neg_dict_elem[0].shape[0]
        H_pos = np.identity(n1) - float(1.0 / n1) * np.ones((n1, n1))
        H_neg = np.identity(n2) - float(1.0 / n2) * np.ones((n2, n2))

        W = (
            pos_sub_mat_kernel.T.dot(H_pos)
            .dot(pos_cl_weigths)
            .dot(H_pos)
            .dot(pos_sub_mat_kernel)
        )
        W += (
            neg_sub_mat_kernel.T.dot(H_neg)
            .dot(neg_cl_weigths)
            .dot(H_neg)
            .dot(neg_sub_mat_kernel)
        )

        #        true_pos_sub_mat_kernel = kernel_sys[0][pos_lab_indices, :]
        #        true_neg_sub_mat_kernel = kernel_sys[0][neg_lab_indices, :]
        #        m1 = true_pos_sub_mat_kernel.shape[0]
        #        m2 = true_neg_sub_mat_kernel.shape[0]
        #        H_true_pos = np.identity(m1) - float(1.0 / m1) * np.ones((m1, m1))
        #        H_true_neg = np.identity(m2) - float(1.0 / m2) * np.ones((m2, m2))
        #        W -= float(1.0 / m1) * true_pos_sub_mat_kernel.T.dot(H_true_pos).dot(true_pos_sub_mat_kernel)
        #        W -= float(1.0 / m2) * true_neg_sub_mat_kernel.T.dot(H_true_neg).dot(true_neg_sub_mat_kernel)

        W *= kernel_sys[0].shape[0]
        W = kernel_sys[2].dot(W).dot(kernel_sys[2])

        svt_eig_vals, svt_eig_vecs = np.linalg.eigh(W)
        return (svt_eig_vals, svt_eig_vecs)

    def __pack_cl_labels(self, y, label_mask, params):
        labels = np.zeros((label_mask.shape[0], params["dim"]))
        for i in range(params["dim"]):
            labels[:, i] = y[label_mask, 0]
        return labels

    def wvar_cl_mode_directions(
        self, sph_quad_eig_sys, y, label_mask, kernel_sys, params
    ):
        labels = self.__pack_cl_labels(y, label_mask, params)
        sph_kb_terms = self.__interpret_cl_constraint(
            kernel_sys, labels, label_mask, params
        )
        mu = self.orth_nu(params, params["dim"], kernel_sys)
        return self.direction_solver.directions(
            kernel_sys, sph_quad_eig_sys, -0.5 * sph_kb_terms[1], params["r"], mu
        )

    def wvar_cl_mode_ho_directions(
        self, sph_quad_eig_sys, y, label_mask, kernel_sys, params
    ):
        labels = self.__pack_cl_labels(y, label_mask, params)
        sph_kb_terms = self.__interpret_cl_constraint(
            kernel_sys, labels, label_mask, params
        )
        return self.direction_solver.hard_orth_directions(
            kernel_sys,
            sph_quad_eig_sys,
            -0.5 * sph_kb_terms[1],
            params["r"],
            params["dim"],
        )

    def cl_mode_directions(self, sph_quad_eig_sys, y, label_mask, kernel_sys, params):
        labels = self.__pack_cl_labels(y, label_mask, params)
        sph_kb_terms = self.__interpret_cl_constraint(
            kernel_sys, labels, label_mask, params
        )
        mu = self.orth_nu(params, params["dim"], kernel_sys)
        return self.direction_solver.directions(
            kernel_sys, sph_quad_eig_sys, -0.5 * sph_kb_terms[1], params["r"], mu
        )

    def cl_mode_ho_directions(
        self, sph_quad_eig_sys, y, label_mask, kernel_sys, params
    ):
        labels = self.__pack_cl_labels(y, label_mask, params)
        sph_kb_terms = self.__interpret_cl_constraint(
            kernel_sys, labels, label_mask, params
        )
        return self.direction_solver.hard_orth_directions(
            kernel_sys,
            sph_quad_eig_sys,
            -0.5 * sph_kb_terms[1],
            params["r"],
            params["dim"],
        )

    def soft_cp_mode_directions(
        self, sph_quad_eig_sys, label_mask, y, kernel_sys, params, const_nu
    ):
        dim = y.shape[1]
        orth_nu = self.orth_nu(params, dim, kernel_sys)
        lin_term = self.__interpret_cp_lin_constraint(
            kernel_sys, y, label_mask, const_nu
        )
        return self.direction_solver.directions(
            kernel_sys, sph_quad_eig_sys, 0.5 * lin_term, params["r"], orth_nu
        )

    def soft_cp_mode_ho_directions(
        self, sph_quad_eig_sys, label_mask, y, kernel_sys, params, const_nu
    ):
        lin_term = self.__interpret_cp_lin_constraint(
            kernel_sys, y, label_mask, const_nu
        )
        return self.direction_solver.hard_orth_directions(
            kernel_sys, sph_quad_eig_sys, 0.5 * lin_term, params["r"]
        )

    def orth_nu(self, params, dim, kernel_sys):
        return -float((kernel_sys[0].shape[0] * params["orth_nu"]) / float(dim))

    def const_nu(self, params, label_mask, kernel_sys):
        n = kernel_sys[0].shape[0]
        x = label_mask.shape[0]
        return float((params["const_nu"] * n) / x)

    def __interpret_cp_lin_constraint(self, kernel_sys, y, label_mask, const_nu):
        lab_sub_K = kernel_sys[0][label_mask, :]
        ellipsoid_lin_term = -2 * const_nu * lab_sub_K.T.dot(y)
        return kernel_sys[2].dot(ellipsoid_lin_term)

    def __interpret_cl_constraint(self, kernel_sys, y, label_mask, params):
        w = self.const_nu(params, label_mask, kernel_sys)
        ellipsoid_lin_term = kernel_sys[0][:, label_mask].dot(y)
        sph_lin_term = w * kernel_sys[2].dot(ellipsoid_lin_term)
        return (None, sph_lin_term)
