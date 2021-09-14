# Build cubic spline for interpolation

# Sources:
# Formulation and implementation are based on:
#   https://www.math.uh.edu/~jingqiu/math4364/spline.pdf or
#   https://web.archive.org/web/20150702075205/https://www.math.uh.edu/~jingqiu/math4364/spline.pdf
#   https://www.uio.no/studier/emner/matnat/math/MAT-INF4130/h17/book2017.pdf chapter 1
#   http://www.thevisualroom.com/tri_diagonal_matrix.html
# Lecture materials
# YouTube Explanation: https://www.youtube.com/watch?v=LaolbjAzZvg
# https://towardsdatascience.com/numerical-interpolation-natural-cubic-spline-52c1157b98ac
# https://people.cs.clemson.edu/~dhouse/courses/405/notes/splines.pdf
# https://timodenk.com/blog/cubic-spline-interpolation/
# Hornbeck, H., 2020. Fast Cubic Spline Interpolation. arXiv preprint arXiv:2001.09253.

PLOT = False
if(PLOT):
    import matplotlib.pyplot as plt
import random
from math import sqrt

class Spline3:
    def __init__(self):
        self.coeff = None

    def fit(self, X, Y):
        self.x = X
        self.y = Y

        delta_x = self._diff(X)
        delta_y = self._diff(Y)

        # Az = d
        # Solving the system for z is based on Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm): Example: https://gist.github.com/cbellei/8ab3ab8551b8dfc8b081c518ccd9ada9

        # Initiate the slices of the matrices
        dim = len(X)
        diagonal_i = [None for _ in range(dim)] # Diagonal Slices from the Triiagonal matrix for i
        diagonal_i1 = [None for _ in range(dim-1)]  # Diagonal Slices from the Triiagonal matrix for i-1
        self.z = [None for _ in range(dim)] # Z coefficients

        # fill diagonals matrices
        diagonal_i[0] = sqrt(2*delta_x[0])
        diagonal_i1[0] = 0.0

        # Solve [L][y] = [B]
        # "The natural spline is defined as setting the second derivative of the first and the last polynomial equal to zero in the interpolation function’s boundary points"
        B0 = 0.0 # natural boundary condition related to the 2nd derivative of the first polynomial = 0
        self.z[0] = B0 / diagonal_i[0]

        for i in range(1, dim-1, 1):
            diagonal_i1[i] = delta_x[i-1] / diagonal_i[i-1]
            diagonal_i[i] = sqrt(2*(delta_x[i-1]+delta_x[i]) - diagonal_i1[i-1] * diagonal_i1[i-1])
            Bi = 6*(delta_y[i]/delta_x[i] - delta_y[i-1]/delta_x[i-1])
            self.z[i] = (Bi - diagonal_i1[i-1]*self.z[i-1])/diagonal_i[i]

        # Last polynomial
        i = dim - 1
        diagonal_i1[i-1] = delta_x[-1] / diagonal_i[i-1]
        diagonal_i[i] = sqrt(2*delta_x[-1] - diagonal_i1[i-1] * diagonal_i1[i-1])
        Bi = 0.0 # # natural boundary condition related to the 2nd derivative of the last polynomial = 0
        self.z[i] = (Bi - diagonal_i1[i-1]*self.z[i-1])/diagonal_i[i]


        # solve [L.T][x] = [y]
        i = dim-1
        self.z[i] = self.z[i] / diagonal_i[i]
        for i in range(dim-2, -1, -1):
            self.z[i] = (self.z[i] - diagonal_i1[i-1]*self.z[i+1])/diagonal_i[i]


        self.coeff = {"x": self.x, "y": self.y, "z": self.z}    # not (x,y,z), Z is not the 3rd dimension in the euclidean space), it is the 2nd derivative of the spline -> coefficients to describe the spline

    def _diff(self, x):
        new_x = [x[i] - x[i-1] for i in range(1, len(x))]
        return new_x

    def compute(self, x, eps=0.00001):
        # Find the nearest neighbors for the interpolated point
        index = 0
        for i in range(len(self.x)):
            if(self.x[i] - x > eps):
                index = i
                break

        x1, x0 = self.x[index], self.x[index-1]   # Neighbours from x-axis
        y1, y0 = self.y[index], self.y[index-1]   # Y-axis correspondeing values to the neighbours
        z1, z0 = self.z[index], self.z[index-1]
        h1 = x1 - x0                             # difference between the x-axis

        # calculate cubic
        y = z0/(6*h1)*(x1-x)**3 + \
            z1/(6*h1)*(x-x0)**3 + \
            (y1/h1 - z1*h1/6)*(x-x0) + \
            (y0/h1 - z0*h1/6)*(x1-x)
        return y

    def get_coeff(self):
        return self.coeff

    def func_out(self, a, b, num_bins=100):
        bin_dim = (b-a)/100
        x = [a+i*bin_dim for i in range(num_bins)]
        bins = [self.compute(i) for i in x]
        return x, bins 
        

n = 0
m = 0
k = 0
x = []
y = []
x_hat = []

# Using input from cli
# n = int(input())
# for i in range(n):
#     x_i = float(input())
#     x.append(x_i)

# m = int(input())
# for i in range(m):
#     y.append([])
#     for j in range(n):
#         y_j_i = float(input()) # y_j^i
#         y[-1].append(y_j_i)

# k = int(input())
# for i in range(k):
#     x_hat_i = float(input())
#     x_hat.append(x_hat_i)

# Using txt file
with open("input.txt", "r") as f:
    lines = f.readlines()

    n = int(lines[0].strip())
    line = lines[1].strip().split(" ")
    for i in range(n):
        x_i = float(line[i])
        x.append(x_i)

    m = int(lines[2].strip())
    for i in range(m):
        line = lines[3+i].strip().split(" ")
        y.append([])
        for j in range(n):
            y_j_i = float(line[j]) # y_j^i
            y[-1].append(y_j_i)

    k = int(lines[3+m].strip())
    line = lines[4+m].strip().split(" ")
    for i in range(k):
        x_hat_i = float(line[i])
        x_hat.append(x_hat_i)

with open("output.txt", "w") as f:
    for i in range(m):
        X = x
        Y = y[i]
        spline3 = Spline3()
        spline3.fit(X, Y)
        outputs = []
        for j in range(k):
            X_hat = x_hat[j]
            out = spline3.compute(X_hat)
            outputs.append(out)
        f.write(" ".join([str(o) for o in outputs]))
        f.write("\n") 
        print(" ".join([str(o) for o in outputs]))
        if(PLOT):
            plt.scatter(X, Y, label="y", marker="*", s=100, c="g")
            plt.scatter(x_hat, outputs, label="y^hat", marker="+", s=100, c="r")

        if(PLOT):
            function_intrep = spline3.func_out(X[0], X[-1])
            plt.plot(*function_intrep, label="Function - Spline3")
            plt.legend()
            plt.show()

postfix = ["st", "nd", "rd"]
for i in range(m):
    error = 0
    out = None
    out_correct = None
    with open("output.txt", "r") as output:
        out = output.readlines()[i].strip().split(" ")
    with open("output_correct.txt", "r") as output:
        out_correct = output.readlines()[i].strip().split(" ")
    errors = [float(out[i]) - float(out_correct[i]) for i in range(k)]
    error = sum(errors)
    mean = error/k
    std = sqrt(sum([(e-mean)**2 for e in errors])/k)
    print(f"Error for {i+1}{postfix[i]} set = {error}\n {errors}\n {mean}+-{std}")