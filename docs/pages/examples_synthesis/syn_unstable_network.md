# Noisy Unstable Network

This example solves a composite optimization problem with three functions
```{math}
\beta^* \in \argmin f(\beta) + \frac{1}{2}\norm{\beta - b_0}_2^2 + \mathbb{I}_{\mathcal{Z}}(\beta)
```

The three functions in the sum are respectively in the classes $S_{m, L}, S_{1, 1}, S_{0, \infty}$. 

The evaluation of $\nabla f$ is inexact $(w = \nabla f(z) + w_p)$, and the transmission  of these noisy gradients occur over an unstable communication channel
```{math}
\begin{align*}
\mat{c}{x^N_{k+1} \hl z^1_k \\ y^1_k} = \mat{cc|cc}{1.2 I & 0 & I & 0\\0 & -0.2 I & 0 & I \hl 0 & I & 0 & 0 \\ I & 0 & 0 & 2 I} \mat{c}{x^N_{k} \hl w^1_k \\ u^1_k}.
\end{align*}
```



Synthesis is performed with $m=1, L=3$. Only the nonsmooth term $\mathbb{I}_{\mathcal{Z}}$ is evaluated implicitly. The overall algorithm satisfy an $\ell_2$-stability specification at minimal rate $\rho$.


 The resulting controller has a worst-case performance of $\rho < 0.9896$. The algorithm is simulated starting at $x_0=0$, for a problem where $f$ is a randomly generated quadratic and $\mathcal{Z}$ is an $L_1$-ball with radius $50$.


<!-- ## Simulation without Noise -->
Figure [1](#fun) plots the function values, and errors over the execution with $w_p = 0$.
:::{figure} _static/unstable_3_state_err_dark.png
:align: center
:class: only-dark
:name: fun
*Figure 1:* Trace of function values and errors
:::

:::{figure} _static/unstable_3_state_err_dark.png
:align: center
:class: only-light
:name: fun
*Figure 1:* Trace of function values and errors
:::


Figure [2](#state) plots the states, oracle input, oracle output over this execution.
:::{figure} _static/unstable_3_state_iter_dark.png
:align: center
:class: only-dark
:name: state
*Figure 2:* Trace of state and oracle evolution
:::

:::{figure} _static/unstable_3_state_iter_dark.png
:align: center
:class: only-light
:name: state
*Figure 2:* Trace of state and oracle evolution (noiseless)
:::



Figure [3](#track) plots the tracking errors of the algorithm based on the Regulator Equation solutions.
:::{figure} _static/unstable_3_tracking_dark.png
:align: center
:class: only-dark
:name: track
*Figure 3:* Trace of tracking errors (noiseless)
:::

:::{figure} _static/unstable_3_tracking_dark.png
:align: center
:class: only-light
:name: track
*Figure 3:* Trace of tracking errors (noiseless)
:::

<!-- ## Simulation with Noise -->
Noise is now injected into the gradient evaluation. The applied performance input $w_p$ is bounded as $\max_{k \in \{0, \ldots, T\}} \norm{w_{p, k}}_2^2 \leq 10$ for all time horizons $T \in \N$.
Figure [4](#fun) plots the persistent tracking and convergence errors with this exogenous noise. 
:::{figure} _static/unstable_noisy_xerr_dark.png
:align: center
:class: only-dark
:name: fun-noise
*Figure 4:* Trace of tracking, optimality, and consensus errors (noisy)
:::

:::{figure} _static/unstable_noisy_xerr_light.png
:align: center
:class: only-light
:name: fun
*Figure 4:* Trace of tracking, optimality, and consensus errors (noisy)
:::


<!-- Figure [5](#state-noise) 
:::{figure} _static/unstable_noisy_err_dark.png
:align: center
:class: only-dark
:name: state-noise
*Figure 5:* Trace of state and oracle evolution (noisy)
:::

:::{figure} _static/unstable_noisy_err_light.png
:align: center
:class: only-light
:name: state-noise
*Figure 5:* Trace of state and oracle evolution (noisy)
::: -->


```{literalinclude} ../../../examples/synthesis/syn_unstable/syn_unstable_cons_robust.m
:linenos: true
:caption: Code for algorithm synthesis with an unstable network
:language: matlab
:lines: 1-34
```