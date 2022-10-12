import numpy as np

def tdglobal(num ,h, esf, Q_g_l, Q_kg):
    #Calculo das tensoes para a lamina k
    tensao_kg = np.zeros(num)
    tensao_kg = tensao_kg.ravel().tolist()
    dk_g = np.dot(Q_g_l, esf.T)
    #Deformações
    deform_g = np.array([[dk_g[0]], [dk_g[1]], [dk_g[2]]])
    deform_g = np.reshape(deform_g, (3,1))
    #Curvatura
    K_g = np.array([[dk_g[3]], [dk_g[4]], [dk_g[5]]])
    K_g = np.reshape(K_g, (3,1))
    def_des = np.zeros(num)
    def_des = def_des.ravel().tolist()
    des = np.zeros(num)
    des = des.ravel().tolist()

    for lamina in range(num):
        des[lamina] = h[lamina]*K_g
        def_des[lamina] = deform_g - des[lamina]
        tensao_kg[lamina] = np.dot(Q_kg[lamina], def_des[lamina])
    return(deform_g, tensao_kg)

    