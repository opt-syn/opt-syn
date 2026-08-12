# Stochastic Gradient Noise

This synthesis example continues {doc}`Analysis with Stochastic Gradient Noise <../examples_analysis/ana_stoch_noise_pgd>`. An algorithm to solve a composite optimization problem $\min_{\beta \in \R^d} f(\beta) + g(\beta)$ must be synthesized with ordering $(\partial g, \nabla f)$. 

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

Synthesis is performed with parameters $m=1, L=10, \Omega = d I$. The controller is synthesized with the settings:
1. $\nabla f$ must be evaluated  explicitly,
2. The worst-case convergence rate must satisfies $\rho < 0.9$,
3. The stochastic sensitivity should be minimized.

The output of Synthesis is a controller with description 

```{math}
\begin{align*}
 \mat{c}{x_{k+1}^1 \\ x_{k+1}^2 \hl  u_k^1 \\ u_k^2} &= \left[\mat{cc|cc}{ 1 & -0.0865 & -0.1815 & -0.0950 \\
0      & 0.0001  & -0.9049 & 0.0950  \hl
1& -0.8144 & -0.8144 & 0       \\
1 & -0.0865 & -0.0865 & 0      } \otimes I_d \right]  \mat{c}{x_{k+1}^1 \\ x_{k+1}^2 \hl   y_k^1 \\ y_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\partial g_(z_k^1) \\  \nabla f(z_k^2)}.
\end{align*}
```

<!-- This controller has a certified convergence rate of $\rho <  0.8561$ and a stochastic sensitivity of $0.3792$.  -->


The controller is used to solve a composite optimization problem in which $f$ is a quadratic and $g$ is the indicator function of an $L_1$ ball. The gradient noise is i.i.d. normally distributed with covariance $100 I$. 


Figure [1](#long-time) plots the squared error  $\norm{z_{p, k}}_2^2$ for 2000 time steps. The empirical mean  (low dotted red line) with value   $150.62$ is computed by averaging the squared error from times 200 to 2000. The mean bound (high dotted gray line) is computed by $d (\text{stdev}^2) \text{sensitivity}^2 =  40 (10^2) (0.3082)^2 =   380.0048$. 

:::{figure} _static/empirical_mean_h2_dark.png
:align: center
:class: only-dark
:name: long-time
*Figure 1* Long-time mean square error bounds
:::

:::{figure} _static/empirical_mean_h2_light.png
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

Figure [3](#noisy-track) plots the tracking errors. After an initial transient, the tracking errors settle to a noise floor.

:::{figure} _static/err_bounded_h2_dark.png
:align: center
:class: only-dark
:name: noisy-track
*Figure 3* Tracking bounds in the noisy algorithm execution
:::

:::{figure} _static/err_bounded_h2_light.png
:align: center
:class: only-light
:name: noisy-track
*Figure 3:* Tracking bounds in the noisy algorithm execution
:::




```{literalinclude} ../../../examples/analysis/syn_quad_h2_single.m
:linenos: true
:caption: Code for Synthesis with stochastic gradient noise
:language: matlab
```

