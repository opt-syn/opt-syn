# Triple Momentum

The Triple Momentum (TMM) algorithm is an explicit first-order accelerated method to optimize a function $f \in S_{m, L}$ {footcite}`van2017fastest`. The worst-case linear convergence rate for the triple momentum scheme is $\rho \leq 1 - \sqrt{\frac{m}{L}}$ .


The Triple Momentum method was extended to composite optimization in {footcite}`upadhyaya2026optimal`. The composite triple momentum algorithm minimizes the sum of functions $f_1 + f_2$  with $f_1 \in S_{m, L}$ and $f_2 \in S_{0, \infty}$, and maintains the worst-case convergence rate $\rho \leq 1 - \sqrt{\frac{m}{L}}$.

Analysis is used to numerically verify the rate  $1 - \sqrt{\frac{m}{L}}$ using IQCs of order `[1, 1]`.  Figure  [1](#tmm-sweep) plots the results of the a parameter sweep. The left plots are the computed upper bounds (solid) v.s. the true rate (dotted lines). The right plots are the error between the computed bounds and the true rate. The top plots report a sweep of $m \in [0.0033
, 1]$ with $L=1$. The bottom plots report a sweep of $L \in [1, 1/0.0033]$ ($L \in [1,303.0303]$). 

:::{figure} _static/tmm_ana_sweep_dark.png
:align: center
:class: only-dark
:name: tmm-sweep
*Figure 2:* Triple Momentum sweeps
:::

:::{figure} _static/tmm_ana_sweep_light.png
:align: center
:class: only-light
:name: tmm-sweep
*Figure 2:* Triple Momentum sweeps
:::

The $m$ sweep is more faithful to the true rate than the $L$ sweep at high values of $m = 1/L$, while the opposite is true for low values of $m$.

```{literalinclude} ../../../examples/analysis/ana_tmm_single.m
:linenos: true
:caption: Code for Triple Momentum analysis at $m=1, L=10$
:language: matlab
```

```{literalinclude} ../../../examples/analysis/ana_tmm_sweep.m
:linenos: true
:caption: Code for Triple Momentum sweep analysis
:language: matlab
:lines: 7-81
```
