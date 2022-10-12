import numpy as np
import math

def matriz(E1, E2, G12, P12, lam_ori_rad, h, num):

    Q_local = np.zeros(num)
    Q_local = Q_local.ravel().tolist()
#---------------------------Matriz de Rigidez----------------------------
    Q11 =  E1**2/(E1-E2*P12**2)
    Q22 = E1*E2/(E1-E2*P12**2)
    Q66 = G12
    Q12 = (P12*E1*E2/(E1-E2*P12**2))
    Q = np.array([[Q11, Q12, 0], [Q12, Q22, 0], [0, 0, Q66]])

#---------------------------Matriz R----------------------------
    R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    R_inv = np.linalg.inv(R)
#---------------------Matriz de Rigidez Transformada---------------------
    for lamina in range(num):
        m = math.cos(lam_ori_rad[lamina])
        n = math.sin(lam_ori_rad[lamina])
        T = np.array([[m**2, n**2, 2*m*n], [n**2, m**2, -2*m*n], [-m*n, m*n, (m**2 - n**2)]])
        T_inv = np.linalg.inv(T)
        Q_local[lamina] = np.linalg.multi_dot([T_inv, Q, R, T, R_inv])

#------------------------------Matriz ABBD------------------------------
    A_h_dif = np.zeros(num)
    B_h_dif = np.zeros(num)
    D_h_dif = np.zeros(num)
    A_prod = np.zeros(num).ravel().tolist()
    B_prod = np.zeros(num).ravel().tolist()
    D_prod = np.zeros(num).ravel().tolist()

    for k in range(num):
        A_h_dif[k] = (h[k+1] - h[k])
        B_h_dif[k] = (h[k+1]**2 - h[k]**2)
        D_h_dif[k] = (h[k+1]**3 - h[k]**3)

        A_prod[k] = A_h_dif[k]*Q_local[k]
        B_prod[k] = B_h_dif[k]*Q_local[k]
        D_prod[k] = D_h_dif[k]*Q_local[k]

    A = sum(A_prod)
    B = sum(B_prod)/2
    D = sum(D_prod)/3
    
    A_ast = np.linalg.inv(A)
    B_ast = - np.dot(A_ast, B)
    C_ast = np.dot(B, A_ast)
    D_ast = D - np.linalg.multi_dot([B, A_ast, B])
    A_l = A_ast + np.linalg.multi_dot([B_ast, np.linalg.inv(D_ast), B_ast.T])
    B_l = np.dot(B_ast, np.linalg.inv(D_ast))
    C_l = - np.dot(D_ast, C_ast)
    D_l = np.linalg.inv(D_ast)
#---------------------------Matriz ABBD Invertida---------------------------
    Q_global = np.array([[A_l[0][0], A_l[0][1], A_l[0][2], B_l[0][0], B_l[0][1], B_l[0][2]], 
                    [A_l[1][0], A_l[1][1], A_l[1][2], B_l[1][0], B_l[1][1], B_l[1][2]], 
                    [A_l[2][0], A_l[2][1], A_l[2][2], B_l[2][0], B_l[2][1], B_l[2][2]],
                    [C_l[0][0], C_l[0][1], C_l[0][2], D_l[0][0], D_l[0][1], D_l[0][2]], 
                    [C_l[1][0], C_l[1][1], C_l[1][2], D_l[1][0], D_l[1][2], D_l[1][2]], 
                    [C_l[2][0], C_l[2][1], C_l[2][2], D_l[2][0], D_l[2][1], D_l[2][2]]])

    return(Q_global,Q_local)
