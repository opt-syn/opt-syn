# Cyclic Coordinate Descent

A coordinate descent algorithm searches over only a subset of variables in each iteration.
{footcite}`wright2015coordinate`.

Letting $i(k)$ be the active coordinate block at iteration $k$, a coordinate descent method to minimize a function $f$ is with stepsize $\gamma>0$ is 
 ```{math}
\begin{align}
    \beta^{i(k)}_{k+1} = \beta^{i(k)}_{k} -\gamma \   [\nabla f(\beta)]^i. 
\end{align}
```

This coordinate descent algorithm is used to solve an unconstrained optimization problem $\min_{\beta \in \R^{90}} f(\beta)$.
The function $f$ is parameterized by scalars $0 < m < W < L$, a constant vector $b$, and a symmetric matrix $Q$ with eigenvalues between $m$ and $W$. The non-quadratic function $f \in S_{m, L}$ is 
```{math}
f(\beta) =\frac{1}{2} \beta^\top Q \beta + b^\top \beta + (L - W)*\log(\1^\top \cosh(\beta)).
```
Coordinate descent with $\gamma = 0.1$ is performed starting from an initial condition $x_0 = 0$. The $c=6$ coordinate blocks are updated cyclically in increasing order.

Figure [1](#coord-6) plots a short trajectory of cyclic coordinate descent.
:::{figure} _static/sim_coord_6_short_dark.png
:align: center
:class: only-dark
:name: coord-6
*Figure 1:* Convergence of coordinate descent algorithm (short time horizon)
:::

:::{figure} _static/sim_coord_6_short_light.png
:align: center
:class: only-light
:name: coord-6
*Figure 1:* Convergence of coordinate descent algorithm (short time horizon)
:::



## Periodic-Orbit Construction

Cyclic coordinate descent schemes can be modeled as {doc}`Periodic-Orbit <../usage/problem_formulation/system/dynamics>` Systems. 
Given a block-size  $c$, we define the permutation matrix 
```{math}
    M_c = \mat{cc}{0 & I_{c-1} \\ 1 & 0}. 
```
The $c$-block cyclic coordinate descent algorithm with $\gamma>0$ is described by

\begin{align*}
\text{Operator}:  & &w_k  = \nabla f(z_k), \\
    \text{Network}: & &  \mat{c}{M^k z_k \\ M^k y_k} &= \left[\mat{cc:cc}{
  0 & 0 & 0 & I_{c-1}  \\
  0 & 0 & 1 & 0 \\
  1 & 0 & 0 & 0 \\
  0 & I_{c-1} & 0 & 0 
    } \otimes I\right]\mat{c}{M^k w_k \\ M^k u_k}, \\
    \text{Gradient Descent}: & & \mat{c}{M^{k+1} x^c_{k+1} \hl M^k u_k} &= \left[\mat{cc|cc}{M I & -\gamma M I \\
    I & 0} \otimes I \right] \mat{c}{M^k x^c_k \hl M^k y_k}.
\end{align*}


Algorithm simulation is performed by using the {class}`opt_system_periodic_orbit` object.

```{literalinclude} ../../../examples/simulation/sim_coord_descent.m
:linenos: true
:caption: Code for coordinate descent with blocksize $c=6$
:language: matlab
```







