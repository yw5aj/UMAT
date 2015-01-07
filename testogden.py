# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 20:36:10 2014

@author: Administrator
"""


import numpy as np
from scipy import linalg


def get_ogden_modulus(f, ogden_param):
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
                    sc_sigma_c += lambda_[a]**(-2) * lambda_[b]**(-2) \
                        * mu * alpha *(-1/3 * lambdabar[a]**alpha 
                        -1/3 * lambdabar[b]**alpha + 1/9 * 
                        (lambdabar**alpha).sum()) * (
                        np.einsum('i, j, k, l', cna.T[a], cna.T[a], cna.T[b],
                                  cna.T[b]))
                # Second term
                for mu, alpha in ogden_param:
                    sc_sigma_c +=  (mu / lambda_[b]**2 * (lambdabar[b]**alpha -
                        (lambdabar**alpha).mean()) - mu / lambda_[a]**2 * (
                        lambdabar[a]**alpha - (lambdabar**alpha).mean())) /\
                        (lambda_[b]**2 - lambda_[a]**2) * (
                    np.einsum('i, j, k, l', cna.T[a], cna.T[b], cna.T[a], 
                              cna.T[b]) +
                    np.einsum('i, j, k, l', cna.T[a], cna.T[b], cna.T[b],
                              cna.T[a]))
            elif a == b:
                # First term
                for mu, alpha in ogden_param:
                    sc_sigma_c += lambda_[a]**(-2) * lambda_[b]**(-2) \
                        * mu * alpha *(1/3 * lambdabar[a]**alpha 
                        + 1/9 * (lambdabar**alpha).sum()) * (
                        np.einsum('i, j, k, l', cna.T[a], cna.T[a], cna.T[b],
                                  cna.T[b]))
    c_sigma_c = 1./det*np.einsum('iI,jJ,kK,lL,IJKL', f,f,f,f,sc_sigma_c)
    return c_sigma_c


def get_modulus(f, ogden_param, d_1, rate='Jaumann', stress='Cauchy'):
    # Get common quantities
    deformation_dict = get_deformation(f)
    det = deformation_dict['det']
    c_sigma_c_iso = get_ogden_modulus(f, ogden_param)
    c_sigma_c_vol = get_vol_modulus(det, d_1, stress='Cauchy', 
                                    rate='Convective')
    c_sigma_c = c_sigma_c_iso + c_sigma_c_vol
    # Return results
    sigma = get_stress(f, ogden_param, d_1, stress='Cauchy')                                  
    c_sigma_j = c_sigma_c + np.einsum('ik, jl', np.eye(3), sigma) + np.einsum(
        'il, jk', sigma, np.eye(3))
    return c_sigma_j


def get_vol_modulus(det, d_1, stress='Cauchy', rate='Jaumann'):
    c_tau_j = 2. * det * (2*det - 1.) / d_1 * np.einsum('ij, kl', np.eye(3),
        np.eye(3))
    c_sigma_j = c_tau_j / det
    sigma = get_vol_stress(det, d_1)
    c_sigma_c = c_sigma_j - np.einsum('ik, jl', np.eye(3), sigma) - \
        np.einsum('jk, il', np.eye(3), sigma)
    if stress == 'Cauchy' and rate == 'Jaumann':
        tangent = c_sigma_j
    elif stress == 'Kirchoff' and rate == 'Jaumann':
        tangent = c_tau_j
    elif stress == 'Cauchy' and rate == 'Convective':
        tangent = c_sigma_c
    return tangent


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
            pk2 += 1. / lambda_i**2 * mu * (lambdabar_i**alpha - 
                np.mean(lambdabar**alpha)) * np.tensordot(cna_i, cna_i, 0)
    # Generate output
    sigma = det**(-1) * f.dot(pk2).dot(f.T)
    if stress == 'Cauchy':
        output = sigma
    elif stress == 'PK2':
        output = pk2
    return output


def get_vol_stress(det, d_1):
    """
    Only Cauchy stress implemented.
    """
    sigma = 2. / d_1 * (det - 1.) * np.eye(3)
    return sigma


def get_deformation(f):
    det = np.linalg.det(f)
    c = f.T.dot(f)
    u = linalg.sqrtm(c)
    lambda_, cna = np.linalg.eig(u)
    fbar = det**(-1/3) * f
    cbar = fbar.T.dot(fbar)
    ubar = linalg.sqrtm(cbar)
    lambdabar, _ = np.linalg.eig(ubar)
    return locals()


if __name__ == '__main__':
    from constants import f
    mu_array, alpha_array = np.array([160e3]), np.array([2.])
    ogden_param = np.c_[mu_array, alpha_array]
    d_1 = .2
    deformation_dict = get_deformation(f)
    det = deformation_dict['det']
    sigma_iso = get_ogden_stress(f, ogden_param, stress='Cauchy')
    sigma_vol = get_vol_stress(det, d_1)
    sigma = sigma_iso + sigma_vol
    c_sigma_j = get_modulus(f, ogden_param, d_1, rate='Jaumann', stress='Cauchy')