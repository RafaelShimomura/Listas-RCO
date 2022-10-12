import numpy as np
import math

def tdlocal(num, lam_ori_rad,tensao_kg, deform_g):
    
    tensao_kl = np.zeros(num)
    tensao_kl = tensao_kl.ravel().tolist()
    deform_kl = np.zeros(num)
    deform_kl = deform_kl.ravel().tolist()

    for lamina in range(num):
        m = math.cos(lam_ori_rad[lamina])
        n = math.sin(lam_ori_rad[lamina])
        T = np.array([[m**2, n**2, 2*m*n], [n**2, m**2, -2*m*n], [-m*n, m*n, (m**2 - n**2)]])

        tensao_kl[lamina] = np.dot(T, tensao_kg[lamina])
        deform_kl[lamina] = np.dot(T, deform_g)
    return(tensao_kl, deform_kl)