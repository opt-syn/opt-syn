<!-- ---

tocdepth: 2

--- -->

# Cocoercive plus Strongly Monotone

This example involves finding the solution to a two-operator inclusion problem, where $F_1$ is $\beta$-cocoercive and $F_1$ is $\mu$-strongly maximal monotone.
The Douglas-Rachford Algorithm with parameters $\lambda = 1, \gamma = 1$ is used to solve this   problem,
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & - I & - I \hl I &- I & 0 \\
 I & -2 I & - I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{F_1(z_k^1) \\  F_2(z_k^2)}.
\end{align*}
```

The explicit minimal rate $\rho$ for this Douglas Rachford algorithm is (Corollary 4.2{footcite}`ryu2020operator`)
\begin{align*}
\rho_{\text{best}} = \begin{cases} 
\abs{1-\frac{\beta}{\beta+1}} & & \beta^2 + \mu \beta + \beta \leq 0 \\
\abs{1-\frac{\1 + \mu\beta}{(\mu+1)(\beta+1)}} & & \mu \beta - \mu - \beta \leq 0 \\
\abs{1 - \frac{\mu}{\mu+1}} & & \mu^2 + \mu \beta + \mu - \beta \leq 0 \\
\frac{1}{2} \frac{\beta + \mu}{\sqrt{\beta \mu(\beta + \mu + 1)}} & & \text{else}. \\
\end{cases}
\end{align*}

Analysis is performed to numerically estimate this minimal rate at order `{1, 1}`. The parameter $\mu$ is fixed to 1, and $\beta$ is swept from  $\beta \in [10^{-2}, 10^2]$. Figure [1](#dr-est-sum) plots the computed and true contraction rates on top, and the error in the estimates on the bottom.


:::{figure} _static/dr_mono_plus_coco_dark.png
:align: center
:class: only-dark
:name: dr-est-sum
:scale: 70%
*Figure 1:* Estimates of the Douglas Rachford Convergence Rate
:::

:::{figure} _static/dr_mono_plus_coco_light.png
:align: center
:class: only-light
:scale: 70%
:name: dr-est-sum
*Figure 1:* Estimates of the Douglas Rachford Convergence Rate
:::


```{literalinclude} ../../../examples/analysis/ana_dr_coco.m
:caption: Code for Douglas-Rachford sweep analysis
:language: matlab
:linenos:  true
:lines: 1-41
```



This example was inspired by the [demonstration](https://autolyap.github.io/examples/douglas_rachford/cocoercive_plus_strongly_monotone/) in AutoLyap.