import numpy as np
from scipy import linalg

F = np.array([[2, 0, 0], [0, 2, 0], [0, 0, .25]])
b = F.dot(F.T)
C = F.T.dot(F)
I = np.eye(3)

lambda2, N = np.linalg.eig(C)

I1 = np.trace(C)
I3 = np.linalg.det(C)

D3 = 2 * lambda2[2] ** 2 - I1 * lambda2[2] + I3 / lambda2[2]

NN3 = lambda2[2] / D3 * (C - (I1 - lambda2[2]) * I +
                         I3 / lambda2[2] * np.linalg.inv(C))

NN3_2 = (C - lambda2[0] * I) / (lambda2[2] - lambda2[0])