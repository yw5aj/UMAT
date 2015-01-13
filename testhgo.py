import numpy as np
from scipy import linalg

from testnh import get_stress_numerical

F = np.eye(3)
F[0, 1] = 2.5e-1*.5
#F[1, 1] = 1.4
#F[0, 0] = 0.5896
#F[2, 2] = 1.2191
#a1 = np.eye(3)[1]
#a2 = np.eye(3)[1]
theta = 15 * np.pi / 180
a1 = np.array((np.cos(theta), np.sin(theta), 0))
a2 = np.array((np.cos(-theta), np.sin(-theta), 0))

C10 = 3.82e3
D = 1e-6
K1 = 996.6e3
K2 = 524.6
Kappa = 0.0

#C10 *= 1e-6
#K1 *= 1e-6
#D *= 1e6

J = np.linalg.det(F)
Fbar = J**(-1./3) * F
C = F.T.dot(F)
Cbar = J**(-2./3) * C
bbar = J**(-2./3) * F.dot(F.T)
Ibar1 = np.trace(Cbar)
U = linalg.sqrtm(C)
R = F.dot(np.linalg.inv(U))


Ibar41 = a1.dot(Cbar.dot(a1))
Ibar42 = a2.dot(Cbar.dot(a2))


Psi_iso = C10 * (Ibar1 -3)
Psi_vol = 1./D * ((J**2-1)/2-np.log(J))

Ebar1 = Kappa * (Ibar1 - 3) + (1 - 3*Kappa) * (Ibar41 - 1)
Ebar2 = Kappa * (Ibar1 - 3) + (1 - 3*Kappa) * (Ibar42 - 1)

Psi_ani = K1 / (2*K2) * (
    np.exp(K2 * (Ebar1 + np.abs(Ebar1))**2/4) + 
    np.exp(K2 * (Ebar2 + np.abs(Ebar2))**2/4) - 2)
Psi = Psi_iso + Psi_ani + Psi_vol
print('Psi: ', Psi)

# Validate stress with numerical method
params_hgo = {
    'model': 'Holzapfel',
    'C10': C10,
    'D': D,
    'K1': K1,
    'K2': K2,
    'Kappa': Kappa,
    'a1': a1,
    'a2': a2}
sigma = get_stress_numerical(F, params_hgo, eps_s=1e-6)
print(sigma)    
def get_sigma_vm(sigma):
    sigma_vm = np.sqrt(1/2*((sigma[0, 0] - sigma[1, 1])**2 + 
        (sigma[1, 1] - sigma[2, 2])**2 + (sigma[2, 2] - sigma[0, 0])**2
        + 6 * (sigma[0, 1]**2 + sigma[1, 2]**2 + sigma[2, 0]**2)))
    return sigma_vm
print(get_sigma_vm(sigma))

# Validate PK2 with absense of the volumetric term
if np.isclose(J, 1, rtol=1e-2): # Following calculatation based on J=1
    II = .5 * (np.einsum('ik, jl -> ijkl', np.eye(3), np.eye(3)) +
               np.einsum('il, jk -> ijkl', np.eye(3), np.eye(3)))
    #psibar41 = (Ibar41-1)*K1*np.exp(K2*(Ibar41-1)**2)
    #psibar42 = (Ibar42-1)*K1*np.exp(K2*(Ibar42-1)**2)
    #S = J**(-2./3) * np.tensordot((II - 1./3 * np.tensordot(np.linalg.inv(C), C, 0
    #    )), 2*C10*np.eye(3), 2) + 2*(psibar41*np.tensordot(a1, a1, 0) + 
    #    psibar42*np.tensordot(a2, a2, 0))
    #sigma = 1/J * F.dot(S).dot(F.T)
    P = II - 1/3*np.tensordot(np.eye(3), np.eye(3), 0)
    tau_tilde = C10 * bbar + 2*K1*Ebar1*np.exp(K2*Ebar1**2) *\
        (Kappa * bbar + (1-3*Kappa)*(np.tensordot(Fbar.dot(a1), Fbar.dot(a1), 0))) +\
        2*K1*Ebar2*np.exp(K2*Ebar2**2) * (Kappa * bbar + 
        (1-3*Kappa)*(np.tensordot(Fbar.dot(a2), Fbar.dot(a2), 0)))
    tau_bar = np.tensordot(P, tau_tilde, 2)
    sigma_bar = tau_bar / J
    sigma_p = np.eye(3) * (J - 1/J) / D
    print(sigma_bar+sigma_p)
    print(get_sigma_vm(sigma_bar))