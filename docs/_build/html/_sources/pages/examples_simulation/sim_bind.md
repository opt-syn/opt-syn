# Repeated Evaluations

This example performs simulation of the 
 the {doc}`Bind <../usage/problem_formulation/system/bind>`. The algorithm with parameters $\lambda, \gamma > 0$ is used to solve a three-operator inclusion problem. 
```{math}
\mat{c}{x_{k+1} \hl z^1_k \\ z^2_k \\ z^3_k \\ z^4_k} &= 
\mat{c|cccc}{I & -2\gamma \lambda I & -\gamma \lambda I& -2\gamma \lambda I  & -\gamma  \lambda I \hl 
I & -\lambda I & 0 & 0 & 0 \\
I & -\lambda I & 0 & 0 & 0 \\
I & -2\lambda I & -\lambda I & -\lambda I & 0 \\
I & -2 \lambda I & -\lambda I & - \lambda I & 0} \mat{c}{x_k \hl w^1_k \\ w^2_k \\ w^3_k \\ w^4_k}, \quad & \mat{c}{w^1_k \\ w^2_k \\ w^3_k \\ w^4_k} &= \mat{c}{F_1 (z^1_k) \\ F_2 (z^2_k)\\ F_3 (z^3_k) \\ F_2(z^4_k)}
```

Operator $F_2$ is evaluated twice per iteration, at instants 2 and 4.
In this simulation, the operators $F_1$ and $F_2$ are gradients of quadratics, and $F_3$ is the subdifferential of an $L_\infty$ norm indicator function. 

Figure [1](#binder) plots an algorithm execution with $\gamma = 0.4, \lambda = 0.25$ starting from a random initial state $x_0$, for a problem with variable  $\beta \in \R^{200}$.

:::{figure} _static/bind_clean_dark.png
:align: center
:class: only-dark
:name: binder
*Figure 1:* Algorithm with bind $[1, 2, 3, 2]$
:::

:::{figure} _static/bind_clean_light.png
:align: center
:class: only-light
:name: binder
*Figure 1:* Algorithm with bind $[1, 2, 3, 2]$
:::

```{literalinclude} ../../../examples/simulation/sim_repeated.m
:linenos: true
:caption: Code for optimization with repeated operator evaluation
:language: matlab
```

