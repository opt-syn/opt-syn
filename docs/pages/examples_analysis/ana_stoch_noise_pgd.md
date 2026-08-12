# Stochastic Gradient Noise

This example performs Analysis of the Projected Gradient Descent (PGD) algorithm to solve the composite optimization problem $\min f(\beta) + g(\beta)$. The gradients of the function $f \in S_{m, L}$ are corrupted by additive stochastic noise. The p.c.c. function $g$ is unaffected by noise. The additive noise $w_{p, k}$ is i.i.d. normally distributed at each time, with matrix covariance $\Omega$. The performance output is the error $z_{p, k} := z^2_k - \beta^*$. The resulting 

The representation of the noise-corrupted PGD algorithm used for Stochastic analysis is 
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2 \\ z_{p, k}} &= \mat{c|cc:c}{I & - \gamma I &  -\gamma I & -\gamma I \hl I &- \gamma I & 0 & 0 \\ I &- \gamma I & 0 & 0 \hdl I &- \gamma I & 0 & 0 }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2 \hdl w_{p, k}}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\partial g_(z_k^1) \\  \nabla f(z_k^2)}.
\end{align*}
```

This representation ensures that the performance input $w$ never enters an oracle evaluation $\partial g$ or $\nabla f$. 

Analysis is performed with paramters $m = 1, L = 10, \Omega = I$. The stepsize parameter $\gamma$ is swept in the range $[0, \frac{2}{L+m}$. The Stochastic Sensitivity and convergence rate are computed for each $\gamma$ using an Analysis program with  orders `[1, 1]` at each oracle.  Figure [1](#pareto) plots the results of two executions of Analysis. The Separate Search poses and solves two separate Analysis. The Joint Search solves a single multi-criteria Analysis program. The conservatism of the joint approach is illustrated by the difference between the curves.

:::{figure} _static/pgd_pareto_h2_stepsize_dark.png
:align: center
:class: only-dark
:name: pareto
*Figure 1* Convergence rate v.s. sensitivity
:::

:::{figure} _static/pgd_pareto_h2_stepsize_light.png
:align: center
:class: only-light
:name: pareto
*Figure 1:* Convergence rate v.s. sensitivity
:::

```{literalinclude} ../../../examples/analysis/ana_pgd_h2.m
:linenos: true
:caption: Code for Analysis with stochastic gradient noise
:language: matlab
```


