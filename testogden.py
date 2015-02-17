# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 20:36:10 2014

@author: Administrator
"""


import numpy as np
from scipy import linalg

delta = np.eye(3)


def get_ogden_modulus_from_ref(f, ogden_param):
    """
    Only Cauchy stress under convective rate implemented.
    """
    # Get common quantities
    deformation_dict = get_deformation(f)
    det = deformation_dict['det']
    lambda_ = deformation_dict['lambda_']
    lambdabar = deformation_dict['lambdabar']
    cna = deformation_dict['cna']
    # Calculate Cauchy sc_sigma_c (for only isochoric part)
    sc_sigma_c = np.zeros((3, 3, 3, 3))
    for a in range(3):
        for b in range(3):
            if a != b:
                # First term
                for mu, alpha in ogden_param:
                    sc_sigma_c += lambda_[a] ** (-2) * lambda_[b] ** (-2) \
                        * mu * alpha * (-1 / 3 * lambdabar[a] ** alpha
                                        - 1 / 3 * lambdabar[b] ** alpha + 1 / 9
                                        * (lambdabar ** alpha).sum()) * (
                        np.einsum('i, j, k, l -> ijkl', cna.T[a], cna.T[a],
                                  cna.T[b], cna.T[b]))
                # Second term
                for mu, alpha in ogden_param:
                    sc_sigma_c += (mu / lambda_[b] ** 2 * (
                        lambdabar[b] ** alpha - (lambdabar ** alpha).mean())
                        - mu / lambda_[a] ** 2 * (
                            lambdabar[a] ** alpha - (lambdabar ** alpha).mean()
                            )) /\
                        (lambda_[b] ** 2 - lambda_[a] ** 2) * (
                        np.einsum('i, j, k, l -> ijkl', cna.T[a], cna.T[b],
                                  cna.T[a], cna.T[b]) +
                        np.einsum('i, j, k, l -> ijkl', cna.T[a], cna.T[b],
                                  cna.T[b], cna.T[a]))
            elif a == b:
                # First term
                for mu, alpha in ogden_param:
                    sc_sigma_c += lambda_[a] ** (-2) * lambda_[b] ** (-2) \
                        * mu * alpha * (1 / 3 * lambdabar[a] ** alpha
                                        + 1 / 9 * (lambdabar ** alpha).sum())\
                        * (np.einsum('i, j, k, l -> ijkl',
                                     cna.T[a], cna.T[a],
                                     cna.T[b], cna.T[b]))
    c_sigma_c = 1. / det * np.einsum('iI,jJ,kK,lL,IJKL -> ijkl',
                                     f, f, f, f, sc_sigma_c)
    return c_sigma_c


def get_ogden_modulus_holzapfel(f, ogden_param):
    """
    Only Cauchy stress under convective rate implemented.
    """
    # Get common quantities
    deformation_dict = get_deformation(f)
    det = deformation_dict['det']
    lambda_ = deformation_dict['lambda_']
    lambdabar = deformation_dict['lambdabar']
    bna = deformation_dict['bna']
    # Calculate Cauchy sc_sigma_c (for only isochoric part)
    c_sigma_c = np.zeros((3, 3, 3, 3))
    for mu, alpha in ogden_param:
        s_iso = lambda_ ** (-2) * mu * (lambdabar ** alpha -
                                        (lambdabar ** alpha).mean())
        for a in range(3):
            for b in range(3):
                if a != b:
                    # First term
                    c_sigma_c += det ** (-1)\
                        * mu * alpha * (-1 / 3 * lambdabar[a] ** alpha
                                        - 1 / 3 * lambdabar[b] ** alpha + 1 / 9 *
                                        (lambdabar ** alpha).sum()) * (
                        np.einsum('i, j, k, l -> ijkl', bna.T[a], bna.T[a],
                                  bna.T[b], bna.T[b]))
                    # Second term
                    c_sigma_c += det ** (-1) * lambda_[a] ** 2 * lambda_[b] ** 2 * (
                        s_iso[b] - s_iso[a]) / (lambda_[b] ** 2 - lambda_[a] ** 2)\
                        * (np.einsum('i, j, k, l -> ijkl', bna.T[a], bna.T[b],
                                     bna.T[a], bna.T[b]) + np.einsum('i, j, k, l -> ijkl',
                                                                     bna.T[a], bna.T[b], bna.T[b], bna.T[a]))
                elif a == b:
                    # First term
                    c_sigma_c += det ** (-1) \
                        * mu * alpha * (1 / 3 * lambdabar[a] ** alpha
                                        + 1 / 9 * (lambdabar ** alpha).sum()) * (
                        np.einsum('i, j, k, l -> ijkl', bna.T[a], bna.T[a],
                                  bna.T[b], bna.T[b]))
    return c_sigma_c


def get_ogden_modulus_paper(f, ogden_param):
    # Get common quantities
    deformation_dict = get_deformation(f)
    det = deformation_dict['det']
    lambda_ = deformation_dict['lambda_']
    lambdabar = deformation_dict['lambdabar']
    bna = deformation_dict['bna']
    i1 = deformation_dict['i1']
    b = deformation_dict['b']
    i3 = det ** 2
    ii = .5 * (np.einsum('ik, jl', delta, delta) +
               np.einsum('il, jk', delta, delta))
    # Calculate beta, gamma and m
    m = [[] for i in range(3)]
    d = [[] for i in range(3)]
    dprime = [[] for i in range(3)]
    dmdg = [[] for i in range(3)]
    beta = [[] for i in range(3)]
    gamma = [[[] for j in range(3)] for i in range(3)]
    ib = .5 * (np.einsum('ac, bd', b, b) + np.einsum('ad, bc', b, b))
    for i in range(3):
        m[i] = np.tensordot(bna.T[i], bna.T[i], 0)
        beta[i] = np.sum(
            [mu * (lambdabar[i] ** alpha - (lambdabar ** alpha).mean())
             for mu, alpha in ogden_param])
        d[i] = 2 * lambda_[i] ** 4 - i1 * \
            lambda_[i] ** 2 + i3 * lambda_[i] ** (-2)
        dprime[i] = 8 * lambda_[i] ** 3 - 2 * i1 * \
            lambda_[i] - 2 * i3 * lambda_[i] ** (-3)
        dmdg[i] = 1 / d[i] * (
            ib - np.tensordot(b, b, 0) + i3 * lambda_[i] ** (-2) * (
                np.tensordot(delta, delta, 0) - ii))\
            + 1 / d[i] * (lambda_[i] ** 2 * (np.tensordot(b, m[i], 0) +
                          np.tensordot(m[i], b, 0)) -
                          1 / 2 * dprime[i] * lambda_[i] *
                          np.tensordot(m[i], m[i], 0))\
            - 1 / d[i] * (i3 * lambda_[i] ** (-2) * (
                np.tensordot(delta, m[i], 0) + np.tensordot(m[i], delta, 0)))
        for j in range(3):
            if not i == j:
                gamma[i][j] = np.sum(
                    [mu * alpha * (-1. / 3 * lambdabar[i] ** alpha -
                     1. / 3 * lambdabar[j] ** alpha
                     + 1. / 9. * (lambdabar ** alpha).sum())
                     for mu, alpha in ogden_param])
            elif i == j:
                gamma[i][j] = np.sum([
                    mu * alpha * (1. / 3 * lambdabar[i] ** alpha +
                                  1. / 9. * (lambdabar ** alpha).sum())
                    for mu, alpha in ogden_param])
    c_sigma_c = np.zeros((3, 3, 3, 3))
    for i in range(3):
        c_sigma_c += 2 / det * beta[i] * dmdg[i]
        for j in range(3):
            c_sigma_c += 1 / det * gamma[i][j] * np.tensordot(m[i], m[j], 0)
    return c_sigma_c


def get_modulus(f, ogden_param, d_1, get_ogden_modulus, rate='Jaumann', stress='Cauchy'):
    # Get common quantities
    deformation_dict = get_deformation(f)
    det = deformation_dict['det']
    c_sigma_c_iso = get_ogden_modulus(f, ogden_param)
    c_sigma_c_vol = get_vol_modulus(det, d_1, stress='Cauchy',
                                    rate='Convective')
    c_sigma_c = c_sigma_c_iso + c_sigma_c_vol
    # Return results
    sigma = get_stress(f, ogden_param, d_1, stress='Cauchy')
    c_sigma_j = c_sigma_c + 0.5 * (np.einsum('ik, jl', delta, sigma) +
                                   np.einsum('jk, il', delta, sigma) + np.einsum('il, jk', delta,
                                                                                 sigma) + np.einsum('jl, ik', delta, sigma))
    return c_sigma_j


def get_vol_modulus(det, d_1, stress='Cauchy', rate='Jaumann'):
    c_tau_j = 2. * det * (2 * det - 1.) / d_1 * np.einsum('ij, kl', delta,
                                                          delta)
    c_sigma_j = c_tau_j / det
    sigma = get_vol_stress(det, d_1)
    c_sigma_c = c_sigma_j - 0.5 * (np.einsum('ik, jl', delta, sigma) +
                                   np.einsum('jk, il', delta, sigma) + np.einsum('il, jk', delta,
                                                                                 sigma) + np.einsum('jl, ik', delta, sigma))
    if stress == 'Cauchy' and rate == 'Jaumann':
        tangent = c_sigma_j
    elif stress == 'Kirchoff' and rate == 'Jaumann':
        tangent = c_tau_j
    elif stress == 'Cauchy' and rate == 'Convective':
        tangent = c_sigma_c
    return tangent


def get_neohookean_modulus(f, nh_param):
    j = np.linalg.det(f)
    b = f.dot(f.T)
    bbar = j ** (-2 / 3) * b
    c10 = nh_param[0]
    ccc_iso = 2 * c10 / j * (1 / 3 * np.trace(bbar) * (
        np.einsum('ik, jl', delta, delta) +
        np.einsum('il, jk', delta, delta) + 2 / 3 *
        np.einsum('ij, kl', delta, delta)) - 2 / 3 * (
        np.einsum('ij, kl', delta, bbar) +
        np.einsum('ij, kl', bbar, delta)))
    return ccc_iso


def get_stress(f, ogden_param, d_1, stress='Cauchy'):
    deformation_dict = get_deformation(f)
    det = deformation_dict['det']
    sigma = get_ogden_stress(f, ogden_param) + get_vol_stress(det, d_1)
    tau = det * sigma
    pk2 = det * np.linalg.inv(f).dot(sigma).dot(np.linalg.inv(f.T))
    if stress == 'Cauchy':
        output = sigma
    elif stress == 'PK2':
        output = pk2
    elif stress == 'Kirchoff':
        output = tau
    return output


def get_ogden_stress(f, ogden_param, stress='Cauchy'):
    """
    Get (Cauchy by default) stress for isochoric part.
    """
    # Get common quantities
    deformation_dict = get_deformation(f)
    det = deformation_dict['det']
    lambda_ = deformation_dict['lambda_']
    lambdabar = deformation_dict['lambdabar']
    cna = deformation_dict['cna']
    # Calculate PK2 stress
    pk2 = np.zeros((3, 3))
    for i, lambdabar_i in enumerate(lambdabar):
        for mu, alpha in ogden_param:
            lambda_i = lambda_[i]
            cna_i = cna.T[i]
            pk2 += 1. / lambda_i ** 2 * mu * (
                lambdabar_i ** alpha - np.mean(lambdabar ** alpha)
                ) * np.tensordot(cna_i, cna_i, 0)
    # Generate output
    sigma = det ** (-1) * f.dot(pk2).dot(f.T)
    if stress == 'Cauchy':
        output = sigma
    elif stress == 'PK2':
        output = pk2
    return output


def get_vol_stress(det, d_1):
    """
    Only Cauchy stress implemented.
    """
    sigma = 2. / d_1 * (det - 1.) * delta
    return sigma


def get_deformation(f):
    det = np.linalg.det(f)
    c = f.T.dot(f)
    i1 = np.trace(c)
    b = f.dot(f.T)
    u = linalg.sqrtm(c)
    r = f.dot(np.linalg.inv(u))
    lambda_, cna = np.linalg.eig(u)
    bna = r.dot(cna)
    fbar = det ** (-1 / 3) * f
    cbar = fbar.T.dot(fbar)
    bbar = fbar.dot(fbar.T)
    ibar1 = np.trace(cbar)
    ubar = linalg.sqrtm(cbar)
    lambdabar = lambda_ * det**(-1/3)
    return locals()


def get_symmetric_part(a):
    map_sym = np.array([[1, 4, 5], [4, 2, 6], [5, 6, 3]]) - 1
    if a.ndim == 2:
        a_sym = np.empty(6)
        for i in range(3):
            for j in range(3):
                a_sym[map_sym[i, j]] = a[i, j]
    elif a.ndim == 4:
        a_sym = np.empty((6, 6))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        a_sym[map_sym[i, j], map_sym[k, l]] = a[i, j, k, l]
    return a_sym


if __name__ == '__main__':
    from constants import f
    f = np.array([[1, 0, 0.45], [0, 1, 0], [0, 0, 1]])
    mu_array, alpha_array = np.array([160e3]), np.array([2.])
    ogden_param = np.c_[mu_array, alpha_array]
    nh_param = np.array((80e3, .2))
    d_1 = .2
    deformation_dict = get_deformation(f)
    det = deformation_dict['det']
    sigma_iso = get_ogden_stress(f, ogden_param, stress='Cauchy')
    sigma_vol = get_vol_stress(det, d_1)
    sigma = sigma_iso + sigma_vol
    c_sigma_j1 = get_modulus(
        f, ogden_param, d_1, get_ogden_modulus_holzapfel, rate='Jaumann',
        stress='Cauchy')
    c_sigma_j2 = get_modulus(
        f, ogden_param, d_1, get_ogden_modulus_paper, rate='Jaumann',
        stress='Cauchy')
    ccc_iso_holzapfel = get_ogden_modulus_holzapfel(f, ogden_param)
    ccc_iso_paper = get_ogden_modulus_paper(f, ogden_param)
    ccc_iso_neohookean = get_neohookean_modulus(f, nh_param)
