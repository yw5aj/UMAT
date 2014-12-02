# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 20:36:10 2014

@author: Administrator
"""


import numpy as np
from scipy import linalg


def get_ogden_stress(dfgrd, ogden_param, output='PK2'):
    # Get common quantities
    deformation_dict = get_deformation(dfgrd)
    det = deformation_dict['det']
    fbar = deformation_dict['fbar']
    lambda_ = deformation_dict['lambda_']
    lambdabar = deformation_dict['lambdabar']
    cna = deformation_dict['cna']
    # Calculate PK2 stress
    pk2 = np.zeros((3, 3))
    for mu, alpha in ogden_param:
        for i, lambdabar_i in enumerate(lambdabar):
            lambda_i = lambda_[i]
            cna_i = cna.T[i]
            pk2 += 1. / lambda_i**2 * mu * (lambdabar_i**alpha - 
                np.mean(lambdabar**alpha)) * np.tensordot(cna_i, cna_i, 0)
    # Generate output
    if output == 'PK2':
        output = pk2
    elif output == 'Cauchy':
        output = det**(-1) * dfgrd.dot(pk2).dot(dfgrd.T)
    return output


def get_vol_stress(det, d_1, output='Cauchy'):
    return 2. / d_1 * (det - 1.) * np.eye(3)


def get_deformation(dfgrd):
    det = np.linalg.det(dfgrd)
    c = dfgrd.T.dot(dfgrd)
    u = linalg.sqrtm(c)
    lambda_, cna = np.linalg.eig(u)
    fbar = det**(-1/3) * dfgrd
    cbar = fbar.T.dot(fbar)
    ubar = linalg.sqrtm(cbar)
    lambdabar, _ = np.linalg.eig(ubar)
    return locals()


if __name__ == '__main__':
    from constants import dfgrd
    mu_array, alpha_array = np.array([2e5]), np.array([2.])
    ogden_param = np.c_[mu_array, alpha_array]
    d_1 = .1
    deformation_dict = get_deformation(dfgrd)
    det = deformation_dict['det']
    sigma_iso = get_ogden_stress(dfgrd, ogden_param, output='Cauchy')
    sigma_vol = get_vol_stress(det, d_1, output='Cauchy')
    print(sigma_iso+sigma_vol)