# Tracking an Oscillator

This example continues the Channel Memory {doc}`simulation <../examples_simulation/sim_channel_memory_symmetric>` demonstration. 

A two-operator problem is solved over a network with channel memory, parameterized by a forgetting factor $\alpha>0$
```{math}
\begin{align*}
z^1_k &= u_k^1 - \alpha z^1_{k-1},  & y_k^1 &= w_k^1 - \alpha y_{k-1}^1, \\
z^2_k &= u_k^2, & y_k^2 &= w_k^2.
\end{align*}
```


An $\alpha$-dependent controller
```{math}
 \mat{c}{x_{k+1}^c \hl u_k^1 \\ u_k^2} =  \mat{c|cc}{I & -\gamma \lambda & \frac{ -\gamma \lambda}{\alpha+1} \hl
 (1+\alpha) I & 0 & 0 \\
 I & -\gamma I & -\frac{ -\gamma }{\alpha+1} I }  \mat{c}{x_{k}^c \hl y_k^1 \\ y_k^2}
```
is used to solve the problem  with values of $\gamma = 0.4, \lambda = 0.2$.
The operator $F_1$ is the subdifferential of a function in $S_{1, 5}$. The operator $F_2$ is the subdifferential of a function in  $S_{0, \infty}$. 

Figure [1](#alpha-sweep) plots Analysis-computed upper-bounds on $\rho$ as the forgetting factor $\alpha$ increases. The same order is used for each operator.   $\rho=2$ is used as an upper bound in bisection: a rate of $\rho<2$ is certified as a worst-case bound.

:::{figure} _static/channel_memory_alpha_dark.png
:align: center
:class: only-dark
:name: alpha-sweep
*Figure 1:* Convergence rate v.s. forgetting factor
:::

:::{figure} _static/channel_memory_alpha_light.png
:align: center
:class: only-light
:name: alpha-sweep
*Figure 1:* Convergence rate v.s. forgetting factor
:::

<!-- The operators $F_1 = \partial f$ and $F_2 = \partial \mathbb{I}_{\text{convex set}}$ of the optimization problem satisfy -->

```{literalinclude} ../../../examples/analysis/ana_delay_memory.m
:caption: Code for Channel Memory sweep analysis
:language: matlab
:linenos:  true
:lines: 1-52
```