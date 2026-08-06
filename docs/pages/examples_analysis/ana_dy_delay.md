# Davis-Yin with Delay

Davis-Yin Splitting is an algorithm to find a zero of the sum of three operators {footcite}`davis2017three`.
The Davis-Yin algorithm with parameters $(\lambda, \gamma)>0$ is described as
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2 \\ z_k^3} &= \mat{c|cc}{I & - \lambda \gamma I & - \lambda \gamma I & - \lambda \gamma I \hl 
 I & - \gamma I & 0 & 0\\
 I & - \gamma I & 0 & 0\\
 I & -2 \gamma I & - I\gamma  & I\gamma}   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2 \\ w_k^3} \in  \mat{c}{F_1(z_k^1) \\  F_2(z_k^2) \\ F_3(z_k^2)}.
\end{align*}
```

This example analyzes the Davis-Yin splitting method under time delays. The operator $F_2$ is evaluated remotely: there is an $h$ step time delay to reach $F_2$, and another $h$ time steps are required to acquire the evaulation of $F_2$. 

All operators are subdifferentials of functions in $S_{m, L}$.
- $F_1$ is p.c.c.
- $F_2 \in S_{1, 1.5}$
- $F_3 \in S_{0, L_3}$ for a paramter $0 < L \leq \infty$. 

The non-delayed Davis-Yin can only achieve linear convergence if $L_3 < \infty$ (Theorem 1, 2 of {footcite}`yi2025convergence`). 

Analysis is performed with orders for Davis-Yin with delays $h \in \{0, \ldots, 4\}$. 
The parameters $\gamma = 1/{L_3}, \lambda = 1$  are used to define the algorithm $L_3$ is swept  in the range $[10^{-2}, 10^{-3}]$. Figure  [1](#dr-delay-opt) plots the outcome of an Analysis run using order  `[2, 0]` at every operator. 
The stability boundary is the line $\rho=1$. The other curves are the convergence rate $\rho$ as a monotonically increasing function of delay: $h=0$ on the bottom in  blue and $h=4$ on top with green. The Davis-Yin algorithm with these parameters is certified as convergent without delay, and is not certified as convergent with delay.

Figure  plots a short trajectory of cyclic coordinate descent.
:::{figure} _static/davis_yin_delay_dark.png
:align: center
:class: only-dark
:name: dr-delay-opt
*Figure 1:* Davis-Yin with $\gamma = 1/{L_3}$ with delays
:::

:::{figure} _static/davis_yin_delay_light.png
:align: center
:class: only-light
:name: dr-delay-opt
*Figure 1:* Davis-Yin with $\gamma = 1/{L_3}$ with delays
:::


Figure  [2](#dr-delay-half) plots the same curves with Davis-Yin parameters $\gamma = 1/(2{L_3}), \lambda = 0.$
:::{figure} _static/davis_yin_delay_half_dark.png
:align: center
:class: only-dark
:name: dr-delay-opt
*Figure 2:* Davis-Yin with $\gamma = 1/(2{L_3})$ with delays
:::

:::{figure} _static/davis_yin_delay_half_light.png
:align: center
:class: only-light
:name: dr-delay-opt
*Figure 2:* Davis-Yin with $\gamma = 1/(2{L_3})$ with delays
:::


```{literalinclude} ../../../examples/analysis/ana_dy_delay.m
:linenos: true
:caption: Code for delayed Davis-Yin analysis at multiple orders, values of $L_3$, and delays $h$
:language: matlab
```


