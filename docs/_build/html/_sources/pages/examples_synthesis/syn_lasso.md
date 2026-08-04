# LASSO

This demonstration applies Synthesis towards LASSO {footcite}`tibshirani1996regression` solvers. The hard LASSO problem is 
\begin{align}
\beta^* \in \argmin_{\norm{\beta}_1 \leq \tau} \frac{1}{2}\norm{E\beta - b}_2^2..
\end{align}

The LASSO problem with $E \in \R^{m \times n}$ is overparameterized if $n > m$, and is  underparameterized if $n < m$. It is assumed that $E$ has full row rank.


In the underparameterized regime, the convex  least-squares function $f(\beta) = \frac{1}{2}\norm{Ex - b}_2^2$  is a member of $S_{0, \sigma_{\max}(E)^2}$. 
In the overparameterized regime the strongly convex $f$ is a member of $S_{\sigma_{\min}(E)^2, \sigma_{\max}(E)^2}$ 

Synthesis of linearly convergent algorithms will take place in the underparameterized regime, and synthesis of ergodically convergent algorithms will occur in the overparameterized regime. These solvers will only implement explicit/forward evaluations of $f$.

## Overparameterized Regime


A LASSO solver is first Synthesized for an overparameterized system with $\sigma_{\max}(E)^2 = 138.90$. The first solver performs explicit evaluations of the least squares cost. Ergodic convergence is approximately confirmed with a constraint violation of $4.68 \times 10^{-12}$. Figure  [1](#lasso-slow) plots a trajectory of this explicit algorithm. 

:::{figure} _static/lasso_under_explicit_long_dark.png
:align: center
:class: only-dark
:name: lasso-slow
*Figure 1:* Explicit overparamterized LASSO solution
:::

:::{figure} _static/lasso_under_explicit_long_dark.png
:align: center
:class: only-light
:name: lasso-slow
*Figure 1:* Explicit overparamterized LASSO solution 
:::

```{literalinclude} ../../../examples/synthesis/syn_lasso_under.m
:caption: Code for Explicit Underparameterized LASSO
:language: matlab
:emphasize-lines: 22
:lines: 1-35
```


## Underparameterized Regime

An explicit algorithm is synthesized for underparameterized setting. The matrix $E \in \R^{130 \times 100}$ has $\sigma_{\min}(E)^2 = 2.1933$ and $\sigma_{\max}(E)^2 = 433.57$. Figure [2](#lasso-rho) plots an execution of the Synthesized algorithm with linear convergence rate $\rho \leq 0.99$. 

:::{figure} _static/lasso_over_explicit_long_dark.png
:align: center
:class: only-dark
:name: lasso-rho
*Figure 2:* Explicit underparamterized LASSO solution 
:::

:::{figure} _static/lasso_over_explicit_long_dark.png
:align: center
:class: only-light
:name: lasso-rho
*Figure 2:* Explicit underparameterized LASSO solution 
:::


```{literalinclude} ../../../examples/synthesis/syn_lasso.m
:caption: Code for Explicit Overparameterized LASSO
:language: matlab
:emphasize-lines: 23-32
:lines:  1-37
```
