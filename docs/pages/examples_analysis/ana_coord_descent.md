# Coordinate Descent Algorithms


This example continues the coordinate-descent {doc}`simulation <../examples_simulation/sim_coord_descent>` demonstration. 

A function $f \in S_{m, L}$ is optimized using a $c$-block cyclic coordinate descent algorithm with gradient steps. 


Analysis is performed for this gradient descent scheme
```{math}
\begin{align}
    \beta^{i(k)}_{k+1} = \beta^{i(k)}_{k} - \gamma \   [\nabla f(\beta)]^i
\end{align}
```
with stepsize $\gamma := \frac{2}{c(m+L)}$. The coordinate block size $c$ is increased from 1 to 6, and the parameter $L$ is swept from 1 to 100 while $m$ is kept to 1. Figure [1](#coord-6) plots the $\rho$ upper bounds computed at order `[1, 0]`.


:::{figure} _static/coord_sweep_ana_dark.png
:align: center
:class: only-dark
:name: coord-6
*Figure 1:* Convergence rate of coordinate descent algorithm v.s. number of blocks
:::

:::{figure} _static/coord_sweep_ana_light.png
:align: center
:class: only-light
:name: coord-6
*Figure 1:* Convergence rate of coordinate descent algorithm v.s. number of blocks.
:::


The specific analysis code at $m=1, L=5, c=4$ is 
```{literalinclude} ../../../examples/analysis/ana_coord_descent.m
:linenos: true
:caption: Code for coordinate descent analysis
:language: matlab
```


:::{seealso}
{doc}`Analysis <../examples_analysis/ana_coord_descent>` of the coordinate-descent scheme.
:::
