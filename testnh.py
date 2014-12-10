# %% Definition section
import numpy as np
from constants import EPS


def perturb(F, i, j, eps=EPS):
    """
    To purturb deformation gradient on its (i, j) component
    """
    e_i, e_j = np.eye(3)[[i, j]]
    del_F = eps/2. * (np.outer(e_i, e_j).dot(F)
        + np.outer(e_j, e_i).dot(F))
    F_hat = F + del_F
    return F_hat


def get_stress_numerical_P(F, params, output='Cauchy'):
    def get_psi(F, params):
        J = np.linalg.det(F)
        Fbar = J**(-1./3.)*F
        Bbar = Fbar.dot(Fbar.T)
        Ibar1 = np.trace(Bbar)        
        if params['model'] == 'Neo-Hookean':
            psi = params['G'] * (Ibar1-3.) + 1./params['D']*(J-1.)**2
        return psi    
    J = np.linalg.det(F)
    P = np.empty((3, 3))
    psi = get_psi(F, params)
    for i in range(3):
        for j in range(3):
            e_i, e_j = np.eye(3)[[i, j]]
            P[i, j] = (get_psi(F + EPS*np.tensordot(e_i, e_j, 0), params) 
                - psi) / EPS
    sigma = 1/J*P.dot(F.T)
    return sigma


def get_stress_numerical_S(F, params, output='Cauchy'):
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
            S[i, j] = (get_psi(C + EPS*(np.tensordot(e_i, e_j, 0) + 
                np.tensordot(e_j, e_i, 0)), params) - psi) / EPS
    sigma = 1/J*F.dot(S).dot(F.T)
    return sigma



def get_stress(F, params, output='Cauchy'):
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


def get_C_CJ_ij(F, params, i, j, eps=EPS):
    J = np.linalg.det(F)
    tau = get_stress(F, params, output='Kirchoff')
    tau_perturb = get_stress(perturb(F, i, j, eps=eps), params, 
                             output='Kirchoff')
    C_CJ_ij = 1. / J / eps * (tau_perturb - tau)
    return C_CJ_ij


def get_C_CJ_numerical(F, params, eps=EPS):
    C_CJ_numerical = np.empty([3, 3, 3, 3])
    for i in range(3):
        for j in range(3):
            C_CJ_numerical[:, :, i, j] = get_C_CJ_ij(F, params, i, j, eps=eps)
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
    params_nh = dict(G=1e5, D=.1, model='Neo-Hookean')
    C_CJ_theoretical = get_C_CJ_theoretical(F, params_nh)
    C_CJ_numerical  = get_C_CJ_numerical(F, params_nh, eps=EPS)
#    print(np.allclose(C_CJ_theoretical, C_CJ_numerical))
    print(get_stress(F, params_nh, output='Cauchy'))
