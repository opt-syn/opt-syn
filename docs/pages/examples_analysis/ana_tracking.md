# Tracking an Oscillator

This example continues the Oscillator {doc}`simulation <../examples_simulation/sim_tracking>` demonstration. 

Figure [1](#track-ana) plots bounds on the convergence rate $\rho$ as $\omega$ is swept in the range $[-\pi, \pi]$. Each Analysis run uses the same orders for each operator. The $\rho=1$ stability boundary is shown by the gray dashed line, convergence is certified if $\rho <1$.

:::{figure} _static/track_circle_2_0_sml_dark.png
:align: center
:class: only-dark
:name: track-ana
*Figure 1:* Convergence rate bound v.s. $\omega$
:::

:::{figure} _static/track_circle_2_0_sml_light.png
:align: center
:class: only-light
:name: track-ana
*Figure 1:* Convergence rate bound v.s. $\omega$
:::


<!-- Figure executes the certified-convergent algorithm at $\omega = \pi/4$. -->
Figure [1](#track-ana) certifies algorithm convergence at $\pi/4 \approx 0.7854$, and does not certify convergence at $\pi/2\approx    1.5708$. Figure [2](#track-ana-conv) plots an algorithm trajectory at the certified-convergent $\omega = \pi/4$. 


:::{figure} _static/track_multi_pi4_dark.png
:align: center
:class: only-dark
:name: track-ana-conv
*Figure 2:* Convergence at $\omega = \pi/4$
:::

:::{figure} _static/track_multi_pi4_dark.png
:align: center
:class: only-light
:name: track-ana-conv
*Figure 2:* Convergence at $\omega = \pi/4$
:::

Figure [3](#track-ana-nonconv) plots a nonconvergent algorithm trajectory at  $\omega = \pi/2$. 


:::{figure} _static/track_multi_pi2_dark.png
:align: center
:class: only-dark
:name: track-ana-nonconv
*Figure 3:* Nonconvergence at $\omega = \pi/2$
:::

:::{figure} _static/track_multi_pi2_light.png
:align: center
:class: only-light
:name: track-ana-nonconv
*Figure 3:* Nonconvergence at $\omega = \pi/2$
:::

```{literalinclude} ../../../examples/analysis/track_analysis_sweep.m
:caption: Code for Oscillator Tracker Analysis at order [2, 0], sweeping over $\omega$
:language: matlab
:linenos:  true
:lines: 1-53
```

