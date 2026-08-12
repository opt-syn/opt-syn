# Stochastic Gradient Noise

This synthesis example continues {doc}`Analysis with Stochastic Gradient Noise <../examples_analysis/ana_stoch_noise_pgd>`. An algorithm to solve a composite optimization problem $\min_{\beta \in \R^d} f(\beta) + g(\beta)$ must be synthesized with ordering $(\partial g, \nabla f)$. The algorithm must evaluate  $\nabla f$ explicitly.



The Operators and Network in the System are 
```{math}
\begin{align*}
 \mat{c}{z_k^1 \\ z_k^2 \\ z_{p, k} \\ y_k^1 \\ y_k^2} &= \mat{cc:c:cc}{ 0 & 0 & 0 & I & 0 \\
0 & 0 & 0 & 0 & I \hdl
0 & 0 & 0 & 0 & I \hdl
I & 0 & 0 & 0 & 0 \\
0 & I & 0 & 0 & 0 \\
 }   \mat{c}{w_k^1 \\ w_k^2 \\ w_{p, k} \\ u_k^1 \\ u_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\partial g_(z_k^1) \\  \nabla f(z_k^2)}.
\end{align*}
```

Synthesis is performed using bisection for $m=1, L=10, \Omega = d I$. The convergence rate $\rho$ is minimized in an outer loop, and the inner loop minimizes the Stochastic Sensitivity. The resulting controller is

```{math}
\begin{align*}
 \mat{c}{x_{k+1}^1 \\ x_{k+1}^2 \hl  u_k^1 \\ u_k^2} &= \left[\mat{cc|cc}{ 1 & -0.1204 & -0.2642 & -0.1438 \\
0      & -0.0011 & -0.9196 & 0.0815  \\
1 & -0.5286 & -0.5286 & 0       \\
1 & -0.1203 & -0.1203 & 0   
 } \otimes I_d \right]  \mat{c}{x_{k+1}^1 \\ x_{k+1}^2 \hl   y_k^1 \\ y_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\partial g_(z_k^1) \\  \nabla f(z_k^2)}.
\end{align*}
```

This controller has a certified convergence rate of $\rho <  0.8561$ and a stochastic sensitivity of $0.3792$. The controller is used to solve a composite optimization problem in which $f$ is a quadratic and $g$ is the indicator function of an $L_1$ ball. The gradient noise is i.i.d. normally distributed and white with standard deviation 10. 


Figure [1](#long-time) plots the squared error  $\norm{z_{p, k}}_2^2$ for 2000 time steps. The empirical mean  (low dotted red line) with value   $150.62$ is computed by averaging the squared error from times 200 to 2000. The mean bound (high dotted gray line) is computed by $d \text{stdev}^2 \text{sensitivity}^2 =  40 (10^2) (0.3792)^2 = 575.02$. 

:::{figure} _static/pgd_pareto_h2_stepsize_dark.png
:align: center
:class: only-dark
:name: long-time
*Figure 1* Long-time mean square error bounds
:::

:::{figure} _static/pgd_pareto_h2_stepsize_light.png
:align: center
:class: only-light
:name: long-time
*Figure 1:* Long-time mean square error bounds
:::

Figure [2](#noisy-exec) plots signals of the  in the first 100 iterations of algorithm execution.

:::{figure} _static/syn_h2_signals_dark.png
:align: center
:class: only-dark
:name: noisy-exec
*Figure 2* Signals in the noisy algorithm execution
:::

:::{figure} _static/syn_h2_signals_light.png
:align: center
:class: only-light
:name: noisy-exec
*Figure 2:* Signals in the noisy algorithm execution
:::



```{literalinclude} ../../../examples/analysis/ana_pgd_h2.m
:linenos: true
:caption: Code for Synthesis with stochastic gradient noise
:language: matlab
```

