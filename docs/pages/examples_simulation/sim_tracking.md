# Tracking an Oscillator

This example involves a time-varying optimization algorithm (from {doc}`Tracking <../usage/problem_formulation/system/tracking>`). The optimal solution $\{\beta^*_k\}$ orbits about a constant point $\beta^*_{\text{center}}$ with a frequency of $\omega = \frac{\pi}{8}$ radians per time step.



The optimization problem is to minimize the sum of four quadratics. Each quadratic $f_{ik}$ shares the same shift,
```{math}
f_{ik}(\beta_k - \beta^*_k) = f_{i0}(\beta_0 - \beta_0) \qquad \qquad \forall k \in \N.
```

The path of the optimal solution obeys known linear dynamics: there exists a state $\eta$ such that 
```{math}
\mat{c}{\eta_{k+1} \hl \beta^*_k } = \mat{ccc}{I & 0 & 0 \\
0 & \cos(\omega) I & -\sin(\omega) I \\
0 & \sin(\omega) I & \cos(\omega) I  \hl
I & I & 0} \eta_k.
```

The algorithm is defined in terms of parameters $b_0 \in \R^3, b_1 \in \R, b_2 \in \R$. An execution of this algorithm with  $\beta \in \R^{4}$ is plotted in Figure [1](#tracker).

:::{figure} _static/track_clean_6_dark.png
:align: center
:class: only-dark
:name: tracker
*Figure 1:* Algorithm with oscillating target
:::

:::{figure} _static/track_clean_6_light.png
:align: center
:class: only-light
:name: tracker
*Figure 1:* Algorithm with oscillating target
:::

The subdifferentials $w$ vectors approach at constant value. The iterates $z$ oscillate according to the moving optimal solution.

```{literalinclude} ../../../examples/simulation/sim_track.m
:linenos: true
:caption: Code for tracking an oscillating target
:language: matlab
```
:::{seealso}
{doc}`Analysis <../examples_analysis/ana_tracking>` of the tracking algorithm.
:::
