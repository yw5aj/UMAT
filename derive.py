from sympy import solve, symbols, diff, simplify, S, init_printing
from sympy.abc import t
from IPython.display import display
import numpy as np


def tensorprod(a, b, out=None):
    """
    Returns the tensor product of two ndarray a and b.
    
    """
    return np.outer(a, b, out=out).reshape(a.shape+b.shape)


if __name__ == '__main__':
    init_printing()
    # Define common variables
    delta = np.full((3, 3), S(1), dtype=object)
    ibar_1, ibar_2 = symbols('Ibar_1, Ibar_2')
    lambdabar_1, lambdabar_2, lambdabar_3 = solve(
        t**3 - ibar_1 * t**2 + ibar_2 * t - 1, t)
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
#    c_sigmaj_iso = np.full((3, 3, 3, 3), S(0))
    c_tauj_iso = \
        2 * (diff(psi_iso, ibar_1) + ibar_1 * diff(psi_iso, ibar_2)) * h_1 +\
        (-2) * diff(psi_iso, ibar_2) * h_2 + 4 * (diff(psi_iso, ibar_1, 2) +
        diff(psi_iso, ibar_2) + 2 * ibar_1 * diff(psi_iso, ibar_1, ibar_2) +
        ibar_1**2 * diff(psi_iso, ibar_2, 2)) * np.outer(bbar, bbar).reshape(
        3, 3, 3, 3) + (-4) * (diff(psi_iso, ibar_1, ibar_2) + ibar_1 * diff(
        psi_iso, ibar_2, 2)) * (bbar.dot(np.outer(bbar, bbar).reshape(3, 3, 3, 
        3)) + (np.outer(bbar, bbar).reshape((3, 3, 3, 3))).dot(bbar)) + 4 * \
        diff(psi_iso, ibar_2, 2) * bbar.dot(np.outer(bbar, bbar).reshape((
        3, 3, 3, 3))).dot(bbar) - 4 / 3 * (diff(psi_iso, ibar_1)  + 2 * ibar_1 
        * diff(psi_iso, ibar_2) + ibar_1 * diff(psi_iso, ibar_1, 2) + 
        (ibar_1**2 + 2 * ibar_2) * diff(psi_iso, ibar_1, ibar_2) + 
        2 * ibar_1 * ibar_2 * diff(psi_iso, ibar_2, 2)) * (
        np.outer(delta, bbar).reshape((3, 3, 3, 3)) + np.outer(bbar, delta
        ).reshape((3, 3, 3, 3))) + 4 / 3 * (2 * diff(psi_iso, ibar_2) 
        + ibar_1 * diff(psi_iso, ibar_1, ibar_2) + 2 * ibar_2
        * diff(psi_iso, ibar_2, 2)) * (np.outer(delta, bbar).reshape(
        3, 3, 3, 3).dot(bbar) + np.outer(bbar.dot(bbar), np.eye(3)).reshape((
        3, 3, 3, 3)))
    c_tauj_vol = det * (det * diff(psi_vol, det, 2) + diff(psi_vol, det)) * \
        np.outer(delta, delta).reshape((3, 3, 3, 3))
    
    # Get c_iso
#    c_iso = 2 / j * ().dot(h_1)
    