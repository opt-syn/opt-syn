# Tracking an Oscillator

This example continues the Oscillator simulation {doc}`../examples_simulation/sim_tracking>` demonstration. 

Figure [1](#track-ana) plots bounds on the convergence rate $\rho$ as $\omega$ is swept in the range $[-\pi, pi]$. Each Analysis run has the same orders for all operator classes. The $\rho=1$ stability boundary is shown by the gray dashed line, convergence is certified if $\rho <1$.

:::{figure} _static/track_circle_2_0_sml_dark.png
:align: center
:class: only-dark
:name: track-ana
:scale: 70%
*Figure 1:* Convergence rate bound v.s. $\omega$
:::

:::{figure} _static/track_circle_2_0_sml_light.png
:align: center
:class: only-light
:scale: 70%
:name: track-ana
*Figure 1:* Convergence rate bound v.s. $\omega$
:::


```{literalinclude} ../../../examples/analysis/track_analysis_sweep.m
:caption: Code for Oscillator Tracker Analysis, sweeping $\omega$
:language: matlab
:linenos:  true
:lines: 1-53
```

