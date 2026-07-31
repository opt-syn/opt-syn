# Time-Varying Delay

This example simulates a Projected Gradient Descent algorithm subject to time-varying delays. We aim to solve the constrained optimization problem

```{math}
\beta^* \in \argmin_{\norm{\beta}_\infty \leq 10} f(\beta),
```
where $f$ is a convex quadratic with eigenvalues between $m=1$ and $L=1.5$. 


Without delays, PGD is described using a parameter $\gamma > 0$ by 
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma I & -\gamma I \hl I &0 & 0 \\
 I & -\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\partial f(z_k^1) \\  \partial I_{\norm{\cdot}_\infty \leq 10}(z_k^2)}
\end{align*}.
```

A time-varying delay $h(k)$ is introduced after  evaluation of $\nabla f$. 
The expression for PGD with time-delays is 
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma I & -\gamma I \hl I &0 & 0 \\
 I & -\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_{k-h(k)}^1 \\ w_k^2}, & \mat{c}{w_{k}^1 \\ w_k^2} \in  \mat{c}{\partial f(z_k^1) \\  \partial I_{\norm{\cdot}_\infty \leq 10}(z_k^2)}
\end{align*}.
```

The delay  $h(k)$ is bounded between 0 and 3 steps at all $h(k)$. It is  temporally restricted according to the logics
:::{list-table}
:header-rows: 1
*  - Type
   - Rule
*  -  Periodic
   -  $h(k)$ increases by 1, and drops to 0 if $h(k)=0$
*  - Snap
   -  $h(k)$ either increases by 1 or drops to 0, but $h(k+1)=0$ always holds if $h(k)=3$.
*  - Contiguous
   -  $h(k)$ increases by 1, stays the same, or decreases by 1,  subject to the constraint $h(k) \in [0, 3]$.
:::

The Snap logic is motivated by communication failures: PGD will use the same stale gradient until a successful transmission is received. Snap and Contiguous are nondeterministic switching logics, while Periodic is deterministic given $h(0)$.

<!-- We model the time-varying delays as a switched system, following the approaches in {footcite}`wen2008switched`, {footcite}`conte2020modeling`.  -->


As an example, a time-varying delay restricted to  $h(k) \in \{0, 1, 2\}$ for each $k \in \N$ can be described using the system matrices
```{math}
\begin{align*}
      \left\{\left(\begin{array}{cc|c}
        0 & 0 &  I \\
        I & 0 &  0 \\ \hline
        0 & 0 &  I
        \end{array}
        \right), \left( \begin{array}{cc|c}
        0 & 0 & I \\
        I & 0 & 0 \\ \hline
        I & 0 & 0
        \end{array}
        \right), \left( \begin{array}{cc|c}
        0 & 0 & I \\
        I & 0 & 0  \\\hline
        0 & I  & 0
        \end{array}
        \right) \right\}. \label{eq:variable_delay_example}
    \end{align*}
```

These delay primitives are used to add delays only on the output of $\partial f$. 

We solve the composite optimization problem using a nominal PGD algorithm with stepsize  $\gamma = 0.05$. All executions are performed starting at $x_0=0$, and at a random initial delay $h(0)$. Figure [1](#delay-per) plots the algorithm trajectory under periodic time-varying delays.

:::{figure} _static/sim_time_var_delay_periodic_dark.png
:align: center
:class: only-dark
:name: delay-per
*Figure 1* Periodic time delays
:::

:::{figure} _static/sim_time_var_delay_periodic_light.png
:align: center
:class: only-light
:name: delay-per
*Figure 1:* Periodic time delays
:::


Figure [1](#delay-snap) plots the algorithm trajectory under periodic time-varying delays.
:::{figure} _static/sim_time_var_delay_snap_dark.png
:align: center
:class: only-dark
:name: delay-snap
*Figure 2* Snap time delays
:::

:::{figure} _static/sim_time_var_delay_snap_light.png
:align: center
:class: only-light
:name: delay-snap
*Figure 2:* Snap time delays
:::


Figure [3](#delay-contiguous) plots the algorithm trajectory under contiguous time-varying delays.
:::{figure} _static/sim_time_var_delay_contiguous_dark.png
:align: center
:class: only-dark
:name: delay-cont
*Figure 3* Contiguous time delays
:::

:::{figure} _static/sim_time_var_delay_contiguous_light.png
:align: center
:class: only-light
:name: delay-cont
*Figure 2:* Contiguous time delays
:::

```{literalinclude} ../../../examples/simulation/sim_time_var_delay.m
:linenos: true
:caption: Code for Projected Gradient Descent with time-varying delays
:language: matlab
```

Introducing a time-varying delay before evaluation of $\partial f$ as $z_{k-h(k)}^1$ instead of $w_{k-h(k)}^1$ fails the Regulator Equation requirement for algorithm convergence.