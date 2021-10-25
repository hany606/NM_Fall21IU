# Assignment 3

Solving System of Linear Algebric equations using Minimal residual method

First of all, simple iteration method: Here, we don't have direct access to A matrix, neither x vector, we have only a multiplaction of A with our vector

x_{k+1} = x_{k} - tau*A@x_{k} + tau*b

tau = 2/(L+l) -> but we don't have the eigen values

Minimal residual method
x_{k+1} = x_{k} - tau_{k}(A@x_{k} - b)

r_{k+1} = r_{k} - tau_{k}(A@r_{k})

tau_{k} = (A@r_{k}, r_{k}) / (A@r_{k}, A@r_{k})


Note: 
* '@' means matrix multiplication, 'b' is 'f' in the slides
* '(a,b)' means that it is equal to a dot product of 'a' and 'b' vectors

## How to compile it?:

```bash
g++ solution.cpp -o solution
```
