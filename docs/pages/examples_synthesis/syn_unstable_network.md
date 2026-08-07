# Unstable Network

This example solves a composite optimization problem with three functions
```{math}
\beta^* \in \argmin f(\beta) + \frac{1}{2}\norm{\beta - b_0}_2^2 + \mathbb{I}_{\mathcal{Z}}(\beta)
```

The three functions in the sum are respectively in the classes $S_{m, L}, S_{1, 1}, S_{0, \infty}$. 

The function $f$ is interfaced over an unstable communication channel
```{math}
\begin{align*}
\mat{c}{x^N_{k} \hl z^1_k \\ y^1_k} = \mat{cc|cc}{1.2 I & 0 & I & 0\\0 & -0.2 I & 0 & I \hl 0 & I & 0 & 0 \\ I & 0 & 0 & 2 I} \mat{c}{x^N_{k} \hl w^1_k \\ u^1_k}.
\end{align*}
```

Synthesis is performed with $m=1, L=3$. Only the nonsmooth term $\mathbb{I}_{\mathcal{Z}}$ is evaluated implicitly. The resulting controller has a worst-case performance of $\rho < 0.9896$. The controller is simulated starting at $x_0=0$, for a problem where $f$ is a randomly generated quadratic and $\mathcal{Z}$ is an $L_1$-ball with radius $50$.


Figure [1](#fun) plots the function values, and errors over the execution.
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
*Figure 2:* Trace of state and oracle evolution
:::



Figure [3](#track) plots the tracking errors of the algorithm, based on the Regulator Equation solutions in `sol.cert.regcl`.
:::{figure} _static/unstable_3_state_iter_dark.png
:align: center
:class: only-dark
:name: track
*Figure 3:* Trace of tracking errors
:::

:::{figure} _static/unstable_3_state_iter_dark.png
:align: center
:class: only-light
:name: track
*Figure 3:* Trace of tracking errors
:::

The computed controller has representation
```{math}
\begin{align*}
A_c := \mat{ccccc}{1.0413  & 0.1080  & -0.0307 & 0.0070  & 0.0498  \\
0.0082  & 0.9481  & -0.0309 & 0.0015  & -0.0023 \\
-0.0225 & -0.1000 & 0.9775  & -0.0039 & -0.0214 \\
2.9000  & 5.7396  & -2.3157 & 1.6939  & 4.3458  \\
0.0433  & 0.1096  & 0.0024  & -0.7070 & -1.1356}, B_c := \mat{ccc}{-0.0172 & -0.0527 & -0.0220 \\
-0.0034 & 0.0382  & 0.0690\\
0.0094  & 0.0305  & 0.0531 \\
-1.2083 & -2.0135 & 0.3022\\
-0.0180 & -0.0170 & -0.0194},
C_c := \mat{ccccc}{1.2000 & 0      & 0      & 0.7217  & 0.9970  \\
1.0000 & 0      & 0      & -0.0145 & -0.0075 \\
1.4757 & 2.7122 & 0.1063 & 0.0854  & 0.3965 },
D_c := \mat{ccc}{
         0         0         0
         0         0         0
   -0.1982   -1.6148   -1.7211}.
\end{align*}
```

```{literalinclude} ../../../examples/simulation/dr_example/syn_unstable_cons.m
:linenos: true
:caption: Code for algorithm synthesis with an unstable network
:language: matlab
```