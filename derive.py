from sympy import solve, symbols, diff, sqrt, simplify, S, init_printing
from sympy.abc import t
from IPython.display import display
import numpy as np


def simplify_all(array):
    """
    Apply `sympy.simplify` to all elements in a ndarray and return a copy.
    """
    shape = array.shape
    new_array = np.empty_like(array.ravel())
    for i, item in enumerate(array.ravel()):
        new_array[i] = simplify(item)
    new_array.reshape(shape)
    return new_array
    


if __name__ == '__main__':
    init_printing()
    # Define common variables
    delta = np.eye(3, dtype=int) * S(1)
    ibar_1, ibar_2 = symbols('Ibar_1, Ibar_2')
    lambda_square = solve(t**3 - ibar_1 * t**2 + ibar_2 * t - 1, t)
    lambdabar_1 = sqrt(lambda_square[0])
    lambdabar_2 = sqrt(lambda_square[1])
    lambdabar_3 = sqrt(lambda_square[2])
    det = symbols('J')
    fbar_11, fbar_12, fbar_13, fbar_21, fbar_22, fbar_23, fbar_31, fbar_32,\
        fbar_33 = symbols('''
        Fbar_11, Fbar_12, Fbar_13, Fbar_21, Fbar_22, Fbar_23, Fbar_31,
        Fbar_32, Fbar_33''')
    fbar = np.array([[fbar_11, fbar_12, fbar_13], [fbar_21, fbar_22, fbar_23], 
                     [fbar_31, fbar_32, fbar_33]])
    bbar = fbar.dot(fbar.T)
    # Define strain energy function
    d_1 = symbols('D_1')
    mu, alpha = symbols('mu, alpha')
    psi_iso = 2 * mu / alpha * (lambdabar_1**alpha + lambdabar_2**alpha
        + lambdabar_3**alpha - 3)
    psi_vol = 1/d_1*(det-1)**2
    # Calculate stress
    dev = np.einsum('ik, jl', np.eye(3, dtype=int), np.eye(3, dtype=int))*S(1)\
        - S(1)/3*np.tensordot(delta, delta, 0)
    p = diff(psi_vol, det)
    sigma_vol = p * delta
    ui1 = diff(psi_iso, ibar_1)
    ui2 = diff(psi_iso, ibar_2)
    sigma_iso = 2/det * np.tensordot(dev, (ui1+ibar_1*ui2)*bbar -
        ui2*bbar.dot(bbar))
    for i, row in enumerate(sigma_iso):
        for j, item in enumerate(row):
            new_item = item.subs([(ibar_1, np.trace(bbar)), (ibar_2, S(1)/2*(
                np.trace(bbar)**2-np.trace(bbar.dot(bbar))))])
            sigma_iso[i, j] = simplify(new_item)
"""    
    # Define H_1
    h_1 = np.full((3, 3, 3, 3), S(0), dtype=object)
    h_2 = np.full((3, 3, 3, 3), S(0), dtype=object)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    h_1[i, j, k, l] = S(1) / 2 * (
                        delta[i, k] * bbar[j, l] + 
                        bbar[i, k] * delta[j, l] +
                        delta[i, l] * bbar[j, k] +
                        bbar[i, l] * delta[j, k])
                    h_2[i, j, k, l] = S(1) / 2 * (
                        delta[i, k] * np.sum([bbar[j, p] * bbar[p, l] 
                            for p in range(3)]) + \
                        delta[j, l] * np.sum([bbar[i, p] * bbar[p, k] 
                            for p in range(3)]) + \
                        delta[i, l] * np.sum([bbar[j, p] * bbar[p, k] 
                            for p in range(3)]) + \
                        delta[j, k] * np.sum([bbar[i, p] * bbar[p, l] 
                            for p in range(3)])) + \
                        bbar[i, k] * bbar[j, l] + bbar[i, l] * bbar[j, k]
    # Get C^sigmaJ_iso
    c_tauj_iso = \
        2 * (diff(psi_iso, ibar_1) + ibar_1 * diff(psi_iso, ibar_2)) * h_1 +\
        (-2) * diff(psi_iso, ibar_2) * h_2 + 4 * (diff(psi_iso, ibar_1, 2) +
        diff(psi_iso, ibar_2) + 2 * ibar_1 * diff(psi_iso, ibar_1, ibar_2) +
        ibar_1**2 * diff(psi_iso, ibar_2, 2)) * np.tensordot(bbar, bbar).reshape(
        3, 3, 3, 3) + (-4) * (diff(psi_iso, ibar_1, ibar_2) + ibar_1 * diff(
        psi_iso, ibar_2, 2)) * (bbar.dot(np.tensordot(bbar, bbar).reshape(3, 3, 3, 
        3)) + (np.tensordot(bbar, bbar, 0)).dot(bbar)) + 4 * \
        diff(psi_iso, ibar_2, 2) * bbar.dot(np.tensordot(bbar, bbar).reshape((
        3, 3, 3, 3))).dot(bbar) - 4 / 3 * (diff(psi_iso, ibar_1)  + 2 * ibar_1 
        * diff(psi_iso, ibar_2) + ibar_1 * diff(psi_iso, ibar_1, 2) + 
        (ibar_1**2 + 2 * ibar_2) * diff(psi_iso, ibar_1, ibar_2) + 
        2 * ibar_1 * ibar_2 * diff(psi_iso, ibar_2, 2)) * (
        np.tensordot(delta, bbar, 0) + np.tensordot(bbar, delta
        , 0)) + 4 / 3 * (2 * diff(psi_iso, ibar_2) 
        + ibar_1 * diff(psi_iso, ibar_1, ibar_2) + 2 * ibar_2
        * diff(psi_iso, ibar_2, 2)) * (np.tensordot(delta, bbar).reshape(
        3, 3, 3, 3).dot(bbar) + np.tensordot(bbar.dot(bbar), np.eye(3)).reshape((
        3, 3, 3, 3)))
    c_tauj_vol = det * (det * diff(psi_vol, det, 2) + diff(psi_vol, det)) * \
       np.tensordot(delta, delta, 0)
"""