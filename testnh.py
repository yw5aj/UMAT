# %% Definition section
import numpy as np
import matplotlib.pyplot as plt
from constants import EPS


def get_stress_numerical(F, params, eps_s=EPS, output='Cauchy'):
    def get_psi(C, params):
        J = np.sqrt(np.linalg.det(C))
        Cbar = J**(-2./3.) * C
        Ibar1 = np.trace(Cbar)   
        if params['model'] == 'Neo-Hookean':
            psi = params['G'] * (Ibar1-3.) + 1./params['D']*(J-1.)**2
        return psi
    J = np.linalg.det(F)
    C = np.dot(F.T, F)
    S = np.empty((3, 3))
    psi = get_psi(C, params)
    for i in range(3):
        for j in range(3):
            e_i, e_j = np.eye(3)[[i, j]]
            psiptb = get_psi(C + eps_s*(np.tensordot(e_i, e_j, 0) + 
                np.tensordot(e_j, e_i, 0)), params)
            S[i, j] = (psiptb - psi) / eps_s
    tau = F.dot(S).dot(F.T)
    sigma = tau / J
    if output == 'Cauchy':
        stress = sigma
    elif output == 'Kirchoff':
        stress = tau
    elif output == 'PK2':
        stress = S
    return stress


def get_stress_theoretical(F, params, output='Cauchy', **keywords):
    """
    Calculates Cauchy stress based on model type.
    """
    J = np.linalg.det(F)
    F_bar = J**(-1./3.) * F
    B_bar = F_bar.dot(F_bar.T)
    if params['model'] == 'Neo-Hookean':
        sigma = 2./J*params['G']*(B_bar - 1./3.*np.trace(B_bar)*np.eye(3)) + \
            2./params['D']*(J-1.)*np.eye(3)
    else:
        raise(NotImplementedError(params['model']))
    # Choose output type
    if output == 'Cauchy':
        stress = sigma
    elif output == 'Kirchoff':
        stress = J * sigma
    return stress


def get_C_CJ_ij(F, params, i, j, eps_c=EPS, eps_s=EPS,
                get_stress=get_stress_numerical):
    def perturb(F, i, j, eps_c=EPS):
        """
        To purturb deformation gradient on its (i, j) component
        """
        e_i, e_j = np.eye(3)[[i, j]]
        del_F = eps_c/2. * (np.outer(e_i, e_j).dot(F)
            + np.outer(e_j, e_i).dot(F))
        F_hat = F + del_F
        return F_hat
    J = np.linalg.det(F)
    tau = get_stress(F, params, output='Kirchoff', eps_s=eps_s)
    tau_perturb = get_stress(perturb(F, i, j, eps_c=eps_c), params, 
                             output='Kirchoff', eps_s=eps_s)
    C_CJ_ij = 1. / J / eps_c * (tau_perturb - tau)
    return C_CJ_ij


def get_C_CJ_numerical(F, params, get_stress=get_stress_numerical, 
                       eps_c=EPS, eps_s=EPS):
    C_CJ_numerical = np.empty([3, 3, 3, 3])
    for i in range(3):
        for j in range(3):
            C_CJ_numerical[:, :, i, j] = get_C_CJ_ij(F, params, i, j, 
                eps_c=eps_c, eps_s=eps_s, get_stress=get_stress)
    return C_CJ_numerical


def get_C_CJ_theoretical(F, params):
    J = np.linalg.det(F)
    F_bar = J**(-1./3.) * F
    B_bar = F_bar.dot(F_bar.T)
    C_CJ_theoretical = np.empty([3, 3, 3, 3])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    if params['model'] == 'Neo-Hookean':
                        C_CJ_theoretical[i, j, k, l] = 2. / J * params['G'] *\
                            (0.5 * (np.eye(3)[i, k] * B_bar[j, l] + 
                            B_bar[i, k] * np.eye(3)[j, l] +
                            np.eye(3)[i, l] * B_bar[j, k] +
                            B_bar[i, l] * np.eye(3)[j, k]) -
                            2./3. * (
                            np.eye(3)[i, j] * B_bar[k, l] + 
                            B_bar[i, j] * np.eye(3)[k, l]) + 
                            2./9. * np.eye(3)[i, j] * np.eye(3)[k, l] *
                            np.trace(B_bar)) + (2. / params['D'] * (2. * J - 1.)
                            * np.eye(3)[i, j] * np.eye(3)[k, l])
    return C_CJ_theoretical


# %% Main code
if __name__ == '__main__':  
    from constants import f as F
    # %% Get related quantities
    params_nh = dict(G=80e3, D=2e-1, model='Neo-Hookean')
    sigma_t = get_stress_theoretical(F, params_nh, output='Cauchy')
    sigma_n = get_stress_numerical(F, params_nh, eps_s=1e-4, output='Cauchy')
    C_CJ_theoretical = get_C_CJ_theoretical(F, params_nh)
    C_CJ_numerical_t  = get_C_CJ_numerical(F, params_nh, get_stress=get_stress_theoretical, eps_c=1e-8, eps_s=1e-8)
    C_CJ_numerical_s  = get_C_CJ_numerical(F, params_nh, get_stress=get_stress_numerical, eps_c=1e-6, eps_s=1e-4)
#    print(np.allclose(C_CJ_theoretical, C_CJ_numerical_s))
    err = np.empty((13, 13))
    for i in range(13):
        for j in range(13):
            eps_c = np.power(10., -1*i)
            eps_s = np.power(10., -1*j)
            C_CJ_theoretical = get_C_CJ_theoretical(F, params_nh)
            C_CJ_numerical_t  = get_C_CJ_numerical(F, params_nh, get_stress=get_stress_theoretical, eps_c=eps_c, eps_s=eps_s)
            C_CJ_numerical_s  = get_C_CJ_numerical(F, params_nh, get_stress=get_stress_numerical, eps_c=eps_c, eps_s=eps_s)
            err[i, j] = np.abs((C_CJ_numerical_s - C_CJ_theoretical)/C_CJ_theoretical).sum()
    plt.imshow(err, vmin=0, vmax=1)
#    fortranccj = np.loadtxt('ccj.txt').reshape((3, 3, 3, 3), order='F')
#    fortrantau = np.loadtxt('tau.txt').reshape((3, 3), order='F')
#    fortranerr = np.abs((fortrandata - C_CJ_theoretical)/C_CJ_theoretical).sum()

