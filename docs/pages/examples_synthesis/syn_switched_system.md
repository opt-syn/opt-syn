# Switched Systems



This is example synthesizes an unconstrained optimization algorithm in which the gradient operator $\nabla f$ is interfaced through a network. The network switches between four subsystems in a ring pattern: at each time step $k$, the mode either stays the same, or increases by 1. An example mode sequence is $(1,1,2,3,3,3,4,4,1,2,3,3,\ldots)$.

The individual subsystems are described by the state-space matrices
```{math}
\begin{align*}
    P_1 &: \mat{cc|cc}{0 & 0 & 1 & 0\\0 & 0.2 & 0 & 1 \hl 0 & 1 & 0 & 0 \\ -1 & 0 & 0 & 0}\otimes I, &   P_2 &: \mat{cc|cc}{0.2 & 0 & 0.25 & 0\\0 & 0.9 & 0 & 1 \hl 0 & 1 & 0 & 0 \\ -0.4 & 0 & -0.5 & 3}\otimes I,  \\
    P_3 &: \mat{cc|cc}{-0.3 & 0 & 0.5 & 0\\0 & -0.5 & 0 & 1 \hl 0 & 1 & 0 & 0 \\ -0.3 & 0 & 0.5 & 1}\otimes I,  &   P_4 &: \mat{cc|cc}{1.2 & 0 & 1 & 0\\0 & -0.2 & 0 & 1 \hl 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 2}\otimes I.
\end{align*}
```
Synthesis is performed for a function $f \in S_{1, 2}$. The returned mode-scheduled controller has a worst-case bound of $\rho \leq 0.9609$ across all switching sequences. 

Figure [1](#fun) plots an execution of the switched optimization algorithm.

:::{figure} _static/ring_uncons_state_dark.png
:align: center
:class: only-dark
:name: fun
*Figure 1:* Trace of algorithm simulation convergence
:::

:::{figure} _static/ring_uncons_state_light.png
:align: center
:class: only-light
:name: fun
*Figure 1:* Trace of algorithm simulation and convergence
:::


Figure [2](#track) plots the  tracking error over this simulation.
:::{figure} _static/ring_uncons_track_dark.png
:align: center
:class: only-dark
:name: track
*Figure 2:* Trace of tracking error
:::

:::{figure} _static/ring_uncons_track_light.png
:align: center
:class: only-light
:name: track
*Figure 2:* Trace of tracking error
:::

```{literalinclude} ../../../examples/synthesis/syn_lse_ring.m
:linenos: true
:caption: Code for ring-switched synthesis at $m=1, L=2$
:language: matlab
:lines: 1-56
```
