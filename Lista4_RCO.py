import numpy as np
import math
from Falha import *
from matriz import matriz
from tdglobal import tdglobal
from tdlocal import tdlocal
def ex5():
    E1 = 77.0e9
    E2 = 75.0e9
    v12 = 0.06
    G12 = 6.50e9
    Xt = 963e6
    Yt = 900e6
    Xc = -856e6
    Yc = -900e6
    S12 = 71e6
    material = np.array([E1, E2, v12, G12, Xt, Yt, Xc, Yc, S12])
    num = 12
    e = 0.29e-3
    Nx = 1000000.00
    Ny = 200000.00
    esf = np.array([[Nx, Ny, 0, 0, 0, 0]])
    teta = np.zeros(num)
    for i in range(num):
        if i % 2:
            teta[i] = 45
        else:
            teta[i] = 0
    
    h_total = num*e
    h = np.zeros(num+1)
    for i in range(np.size(h)):
        h[i] = -h_total/2 + e*i

    lam_ori = np.array([teta[0], teta[1], teta[2], teta[3], teta[4], teta[5], teta[6], teta[7], teta[8], teta[9], teta[10], teta[11]])
    teta_rad = teta*math.pi/180
    lam_ori_rad = np.array([teta_rad[0], teta_rad[1], teta_rad[2], teta_rad[3], teta_rad[4], teta_rad[5], teta_rad[6], teta_rad[7], teta_rad[8], teta_rad[9], teta_rad[10], teta_rad[11]])

    [Q_global_5, Q_local_5] = matriz(E1, E2, G12, v12, lam_ori_rad, h, num)
    [deform_g_5, tensao_lam_g_5] = tdglobal(num ,h, esf, Q_global_5, Q_local_5)
    [tensao_lam_l_5, deform_lam_l_5] = tdlocal(num, lam_ori_rad,tensao_lam_g_5, deform_g_5)

    MS_TH = TsaiHill(tensao_lam_g_5, material)
    MS_mais, MS_menos = TsaiWu(tensao_lam_g_5, material)

    return(MS_TH, MS_mais, MS_menos, MS_max)
