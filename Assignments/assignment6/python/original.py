# Source: https://scipython.com/book2/chapter-8-scipy/examples/solving-a-system-of-stiff-odes/
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

k1, k2, k3, k4, k5 = 1.28, 2400000, 33.6, 2400, 1
X0, Y0, Z0, A0, B0, P0 = 0.004, 0.0001, 0.08, 10, 10, 0


def deriv(t, y):
    """ODEs for Robertson's chemical reaction system."""
    # x, y, z = y
    X,Y,Z,A,B,P = y
    xdot = k1*A*Y - k2*X*Y + k3*A*X - 2*k4*X*X #-0.04 * x + 1.e4 * y * z
    ydot = -k1*A*Y - k2*X*Y + k5*B*Z #0.04 * x - 1.e4 * y * z - 3.e7 * y**2
    zdot = 2*k3*A*X - k5*B*Z #3.e7 * y**2
    pdot = k1*A*Y + 2*k2*X*Y + k4*X*X
    adot = -k1*A*Y - k3*A*X + k4*X*X
    bdot = -k5*B*Z
    return xdot, ydot, zdot, adot, bdot, pdot

# Initial and final times.
t0, tf = 0, 500
# Initial conditions: [X] = 1; [Y] = [Z] = 0.
y0 = [X0, Y0, Z0, A0, B0, P0]
# Solve, using a method resilient to stiff ODEs.
soln = solve_ivp(deriv, (t0, tf), y0, method='Radau')
print(soln.nfev, 'evaluations required.')


print(dict(zip(["X","Y","Z","A","B","P"],soln.y[:, -1])))
# Plot the concentrations as a function of time. Scale [Y] by 10**YFAC
# so its variation is visible on the same axis used for [X] and [Z].
plt.plot(soln.t, soln.y[0], label='[X]')
plt.plot(soln.t, soln.y[1], label='[Y]')
plt.plot(soln.t, soln.y[2], label='[Z]')
plt.plot(soln.t, soln.y[3], label='[A]')
plt.plot(soln.t, soln.y[4], label='[B]')
plt.plot(soln.t, soln.y[5], label='[P]')
plt.xlabel('time /s')
plt.ylabel('concentration /arb. units')
plt.legend()
plt.show()