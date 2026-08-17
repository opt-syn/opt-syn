# Noisy Douglas Rachford

The Douglas-Rachford algorithm is a procedure for solving a two-operator inclusion problem {footcite}`douglas1956numerical`.
It is characterized by parameters $\gamma, \lambda > 0$, and can be described as the interconnection
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma \lambda I & -\gamma \lambda I \hl I &-\gamma I & 0 \\
 I & -2\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{F_1(z_k^1) \\  F_2(z_k^2)}.
\end{align*}
```

The Douglas-Rachford algorithm is used in this example to solve a composite optimization problem
```{math}
\beta^* \in \argmin_{\norm{\beta}_\infty \leq 10} f(\beta)
```

The necessary optimality condition for this problem is posed using the operators $F_1 = \partial f$, and $F_2 = \partial \mathbb{I}_{\norm{\cdot}_\infty \leq 10}$ 
```{math}
0 \in \partial  f(\beta^*) + \partial \mathbb{I}_{\norm{\cdot}_\infty}(\beta^*).
```

## No Noise
The Douglas-Rachford scheme with parameters $\gamma = 0.4, \lambda = 1$ is executed for a problem where $f$ is a convex quadratic (eigenvalue bounds $m = 1, L = 10$). Figure [1](#dr-clean) plots a trajectory of Douglas-Rachford starting from $x_0 = 0$.

:::{figure} _static/dr_clean_dark.png
:align: center
:class: only-dark
:name: dr-clean
*Figure 1:* Convergence without noise
:::

:::{figure} _static/dr_clean_light.png
:align: center
:class: only-light
:name: dr-clean
*Figure 1:* Convergence without noise
:::


## With Noise
Noise is then added to the Douglas-Rachford execution. The performance input $w_p$ introduces additive noise
at the output of the subgradient evaluations. The performance output $z_p$ is the consensus error $\pm \frac{1}{2}(z_1 - z_2)$. 
The {doc}`System <../usage/problem_formulation/system/index_system>` representing Douglas-Rachford with this noise structure is
```{math}
\begin{align*}
 \text{Operator} & & \mat{c}{w_k^1 \\ w_k^2} &\in  \mat{c}{F_1(z_k^1) \\  F_2(z_k^2)}, \\
\text{Network} & & \mat{c}{z_k^1 \\ z_k^2 \hdl z_{p, k}^1 \\ z_{p, k}^2 \hdl y_{k}^1\\y_k^2} &= \mat{cc:cc:cc}{0 & 0 & 0 & 0 & I & 0 \\
0 & 0 & 0 & 0 & 0 & I  \hdl
0 & 0 & 0 & 0 & \frac{1}{2} I & -\frac{1}{2} I \\
0 & 0 & 0 & 0 & -\frac{1}{2} I & \frac{1}{2} I \hdl
0 & 0 & I & 0 & I & 0 \\
0 & 0 & 0 & I & 0 & I} \mat{c}{w_k^1 \\ w_k^2 \hdl w_{p, k}^1 \\ w_{p, k}^2 \hdl u_{k}^1 \\ u_k^2}, \\
 \text{Douglas-Rachford} & & \mat{c}{x^c_{k+1} \hl u_k^1 \\ u_k^2} &= \mat{c|cc}{I & -\gamma \lambda I & -\gamma \lambda I \hl I &-\gamma I & 0 \\
 I & -2\gamma I & -\gamma I  }   \mat{c}{x^c_{k} \hl y_k^1 \\ y_k^2}.
\end{align*}
```

Figure [2](#dr-noisy) plots a trace of a trajectory starting at $x_0=0$, in which the performance input $w_p$ is randomly sampled subject to the bound $\norm{w_{p,k}}_2 \leq 10, \forall k \in \N$.
:::{figure} _static/dr_noisy_dark.png
:align: center
:class: only-dark
:name: dr-noisy
*Figure 2:* Response under bounded noise
:::

:::{figure} _static/dr_noisy_light.png
:align: center
:class: only-light
:name: dr-noisy
*Figure 2:* Response under bounded noise
:::


```{literalinclude} ../../../examples/simulation/dr_example/sim_dr_opt_noise.m
:linenos: true
:caption: Code for Douglas-Rachford with noise corruption
:language: matlab
```

