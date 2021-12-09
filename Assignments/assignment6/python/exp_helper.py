print("A_eig <<", end="")
for i in range(6):
    for j in range(6):
        print(f" A_mat[{i*6+j}],", end="")    


# mat = [["-k[1]*u[1]+k[2]*u[3]-2*k[3]*u[0]", "k[0]*u[3]", "0", "0", "0", "0"], 
#               ["0", "-k[1]*u[0]-k[0]*u[3]", "k[4]*u[4]", "0", "0", "0"],
#               ["0", "0", "-k[4]*u[4]", "2*k[2]*u[0]", "0", "0"],
#               ["k[3]*u[0]", "0", "0", "-k[0]*u[1]-k[2]*u[0]", "0", "0"],
#               ["0", "0", "0", "0", "-k[4]*u[2]", " 0"],
#               ["k[3]*u[0]", "2*k[1]*u[0]", "0", "k[0]*u[1]", "0", "0"]]

# for i in range(6):
#     for j in range(6):
#         print(f"A[{i*6+j}] = {mat[i][j]};")    
#     print("\n")  



# import sympy as sp
# # from math import sin, cos
# import numpy as np

# def print_mat(mat):
#     for i in range(6):
#         for j in range(6):
#             print(mat[i,j],end=",\t")
#         print("\n")
#     print("\n\n\n\n-------------\n\n\n\n")

# mat = sp.Matrix
# x, y, z, a, b, p = sp.symbols("x y z a b p", real=True)
# k1, k2, k3, k4, k5, k6 = sp.symbols("k1 k2 k3 k4 k5 k6", real=True)

# A = mat(\
#     np.array([[-k2*y+k3*a-2*k4*x, k1*a, 0, 0, 0, 0],
#               [0,                -k2*x-k1*a, k5*b,    0, 0, 0],
#               [0,                   0,       -k5*b, 2*k3*x, 0, 0],
#               [k4*x,                   0,         0,  -k1*y-k3*x, 0, 0],
#               [0, 0, 0, 0, -k5*z, 0],
#               [k4*x, 2*k2*x, 0, k1*y, 0, 0]])
#         )


# lam = sp.symbols('lambda')
# cp = sp.det(A - lam * sp.eye(6))
# eigs = sp.roots(sp.Poly(cp, lam))
# # eigs = A.eigenvals()
# print(eigs)

# print_mat(A)
# print("\n\n\n\n-------------\n\n\n\n")

# print(sp.printing.latex(A))

# # q1, q2, q3, q4, q5, q6 = 0,0,0,0,0,0
# # l1, l2, l3, l4, l5, l6 = 0,0,0,0,0,0

# # For position vector: in order to get q1, q2, q3
# # x = l2*cos(q1) + l3*cos(q1)*cos(q2) + l4*(-sin(q2)*sin(q3)*cos(q1) + cos(q1)*cos(q2)*cos(q3)) + l5*(-sin(q2)*sin(q3)*cos(q1) + cos(q1)*cos(q2)*cos(q3))
# # y = l2*sin(q1) + l3*sin(q1)*cos(q2) + l4*(-sin(q1)*sin(q2)*sin(q3) + sin(q1)*cos(q2)*cos(q3)) + l5*(-sin(q1)*sin(q2)*sin(q3) + sin(q1)*cos(q2)*cos(q3))
# # z = -l3*sin(q2) + l4*(-sin(q2)*cos(q3) - sin(q3)*cos(q2)) + l5*(-sin(q2)*cos(q3) - sin(q3)*cos(q2))

# # Full matrix to substitute in it q1, q2, q3 after getting them using the last column
# # In order to be able to get its inverse and get T_{456}_{bits} = T_{123}^{-1} T_o = T_{456}
# # T_123 = [
# #             [-sin(q2)*sin(q3)*cos(q1) + cos(q1)*cos(q2)*cos(q3),	-sin(q1),	sin(q2)*cos(q1)*cos(q3) + sin(q3)*cos(q1)*cos(q2),	l2*cos(q1) + l3*cos(q1)*cos(q2) + l4*(-sin(q2)*sin(q3)*cos(q1) + cos(q1)*cos(q2)*cos(q3)) + l5*(-sin(q2)*sin(q3)*cos(q1) + cos(q1)*cos(q2)*cos(q3)),]
# #             [-sin(q1)*sin(q2)*sin(q3) + sin(q1)*cos(q2)*cos(q3),	cos(q1),	sin(q1)*sin(q2)*cos(q3) + sin(q1)*sin(q3)*cos(q2),	l2*sin(q1) + l3*sin(q1)*cos(q2) + l4*(-sin(q1)*sin(q2)*sin(q3) + sin(q1)*cos(q2)*cos(q3)) + l5*(-sin(q1)*sin(q2)*sin(q3) + sin(q1)*cos(q2)*cos(q3)),]
# #             [-sin(q2)*cos(q3) - sin(q3)*cos(q2),	0,	-sin(q2)*sin(q3) + cos(q2)*cos(q3),	-l3*sin(q2) + l4*(-sin(q2)*cos(q3) - sin(q3)*cos(q2)) + l5*(-sin(q2)*cos(q3) - sin(q3)*cos(q2)),]
# #             [0,	0,	0,	1,]
# #         ]

# # For rotation matrix: in order to get q4, q5, q6
# # cos(q5),	sin(q5)*sin(q6),	sin(q5)*cos(q6)
# # sin(q4)*sin(q5),	-sin(q4)*sin(q6)*cos(q5) + cos(q4)*cos(q6),	-sin(q4)*cos(q5)*cos(q6) - sin(q6)*cos(q4)	
# # -sin(q5)*cos(q4),	sin(q4)*cos(q6) + sin(q6)*cos(q4)*cos(q5),	-sin(q4)*sin(q6) + cos(q4)*cos(q5)*cos(q6)