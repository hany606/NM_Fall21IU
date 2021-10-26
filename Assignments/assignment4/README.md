# Assignment 4

Solving system of non-linear equations

Newton's method
Source: 
* https://personal.math.ubc.ca/~pwalls/math-python/roots-optimization/newton/
* https://www.lakeheadu.ca/sites/default/files/uploads/77/docs/RemaniFinal.pdf


## Method:
* F = [f_1; f_2; f_3; f_4; f_5, ...., f_n].T
* \phi = \sum_{i=0}^{n} r_i(x,y,z,t)       (Lecture notations)
* r_i(x,y,z,t) = f_i
* f_i = (x-x_i)^2 + (y-y_i)^2 + (z-z_i)^2 - (t-t_i)^2 = 0 
* Convert the problem to optimization min_{x,y,z,t} (\phi) -> this means we will satisfy the equations and get the values of x,y,z,t that will make the summation at least as possible (satisfy the equality)
* Solve it using Gradient based optimization method for unconstrained optimization, we might add a constraint on time to be positive but we can change it back to unconstrained by adding this constraint as extra cost (Lagrange multiplier)
* Moreover, we can solve using simple steepest gradient descent or Newton's method:
 1. Get the jacobian matrix (1st order derivative) (nx4) such that the row: [\partial f_i / \partial x, \partial f_i / \partial y, \partial f_i / \partial z, \partial f_i / \partial t] \forall i \in {0..n}
 2. Get the inverse of the jacobian using pseudo inverse as it is not squared matrix
 3. Multiple the inverse of the jacobian with F   (Update step)
 4. Subtract the update step vector from the previous iteration result
* And repeat these steps for each iteration till convergence with error margin for each new incoming data


Notes:
* Gradient matrix (Jacobian)  J: nx4 -> s.t. n is the number of the satalities that we recieved their data
* for n = 4 -> J: 4x4, we will just calculate the inverse
* otherwise; n > 4 -> J: nx4, we will calculate the pseudo inverse (J^T@J)^{-1}@J^T
*   Thus, we will need inverse of (J^T@J)
