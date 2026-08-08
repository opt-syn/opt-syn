# Proximal Heavy Ball

The Proximal Heavy Ball algorithm is a scheme to solve composite optimization problems $\min_\beta f(\beta) + g(\beta)$ (Sec. 2.5.3. {footcite}`upadhyaya2025automated`). 


solves composite optimization problems with $f \in S_{m, L}$ and $g \in S_{0, \infty}$. The description of the Proximal Heavy Ball algorithm is with parameters $\gamma, \lambda > 0$ is
```{math}
\begin{align*}
 \mat{c}{x_{k+1}^1 \\ x_{k+1}^2 \hl z_k^1 \\ z_k^2} &= \mat{cc|cc}{ (1+\lambda) I & -\lambda I & -\gamma I & -\gamma I \\
  I & 0 & 0 & 0 \hl
 I & 0 & 0 & 0 \\
 (1+\lambda) I & -\lambda I & -\gamma I & -\gamma I }   \mat{c}{x_{k}^1 \\ x_k^2 \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\nabla f(z_k^1) \\  \partial g(z_k^2)}.
\end{align*}
```

Analysis is performed to upper-bound the convergence worst-case rate $\rho$ for parameter choices $\gamma = \frac{2}{m+L}$ and $\lambda = 0.65$. 
The operators $\nabla f$ and $\partial g$ are Analyzed each with orders `[3, 0]` or `[3, 1]`.
The  parameter $L$ swept in the range $[1, 500]$ with fixed $m=1$. $f$ is restricted to general functions in $S_{m, L}$, or to quadratics in $S_{m, L}$.


Figure [1](#phb-quad) plots the computed $\rho$ upper bounds at each order and operator class. The worst-case convergence rate $\rho$ is reduced by restricting to quadratic $f$ as compared to general $f \in S_{m, L}$.

:::{figure} _static/ana_phb_quad_dark.png
:align: center
:class: only-dark
:name: phb-quad
*Figure 1:* Estimates of the Proximal Heavy Ball convergence rate
:::

:::{figure} _static/ana_phb_quad_light.png
:align: center
:class: only-light
:scale: 70%
:name: phb-quad
*Figure 1:* Estimates of the Proximal Heavy Ball convergence rate
:::


```{literalinclude} ../../../examples/analysis/ana_phb_sweep.m
:caption: Code for Proximal Heavy Ball sweep analysis
:language: matlab
:linenos:  true
```


