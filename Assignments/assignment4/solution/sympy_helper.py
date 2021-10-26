# Work too slow to generate the analytical values for the gradients
import sympy as sym

N = 30

def print_symbols():
    symbols = [f"x{i}" for i in range(N)]
    
    print(f"Symbols string:")
    print(' '.join(symbols))

    print(f"vars string:")
    print(', '.join(symbols))

def define_matrices(n):
    prefix = "sym.Matrix(["
    postfix = "])"

    rows = []

    for i in range(n):
        rows.append(f"[2*(x-x{i}), 2*(y-y{i}), 2*(z-z{i}), 2*(t-t{i})]")

    output = prefix + ",\n".join(rows) + postfix
    return output
# print_symbols()

x, y, z, t = sym.symbols("x y z t")

x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29 = sym.symbols("x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29")
y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, y21, y22, y23, y24, y25, y26, y27, y28, y29 = sym.symbols("y0 y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16 y17 y18 y19 y20 y21 y22 y23 y24 y25 y26 y27 y28 y29")
z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24, z25, z26, z27, z28, z29 = sym.symbols("z0 z1 z2 z3 z4 z5 z6 z7 z8 z9 z10 z11 z12 z13 z14 z15 z16 z17 z18 z19 z20 z21 z22 z23 z24 z25 z26 z27 z28 z29")
t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29 = sym.symbols("t0 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 t13 t14 t15 t16 t17 t18 t19 t20 t21 t22 t23 t24 t25 t26 t27 t28 t29")

# sx, sy, rho = sym.symbols('sigma_x sigma_y rho')



matrix = define_matrices(n=2)
print(matrix)

matrix = eval(matrix)
print(matrix)
print(matrix.pinv())