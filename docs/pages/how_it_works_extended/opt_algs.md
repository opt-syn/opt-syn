## Optimization/Inclusion Problems


An optimization problem tries to find a point $\beta^*$ satisfying 
```{math}
\beta^* & \in \argmin_\beta \sum_{i=1}^s f_i(\beta).
```

An inclusion problem is a more general concept than an optimization problem. Inclusion problems try to find a point $\beta^*$ in the zero in the sum of operators.
```{math}
\begin{align*}
 0 \in \sum_{i=1}^s F_i(\beta^*).
\end{align*}
```

A solution to a zero-inclusion problem is a pair $(\beta^*, w^*)$ with 

```{math}
\begin{align*}
 0 \in \sum_{i=1}^s w^{*i}, \qquad w^{*i} \in F_i(\beta^*).
\end{align*}
```

Zero-inclusion problems include solution concepts such as optimization problems and variational inequalities. In particular, if each $F_i$ is the subdifferential of a proper function $f_i$, then 
```{math}
\begin{align*}
 \beta^* & \in \argmin_\beta \sum_{i=1}^s f_i(\beta) &  \text{implies} & & 0 & \in \sum_{i=1}^s  F_i(\beta^*).
\end{align*}
```
The zero inclusion problem in $\partial f$ is an optimality principle for the optimization problem.




## Optimization Algorithms



An optimization algorithm is a procedure that generates a sequence of iterates $(w_k, z_k)_{k \in \N}$ satisfying $w^i_k \in F^i(z^i_k)$.


 Many common optimization algorithms can be expressed as the interconnection of operators and linear systems. As an example, the gradient descent/forward-step method with stepsize $\gamma > 0$ may be represented by 
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k} &= \mat{c|cc}{I & -\gamma \lambda I \hl I & 0 }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{F_1(z_k^1)},
\end{align*}
```
 and the   Douglas-Rachford algorithm {footcite}`douglas1956numerical`   with parameters $\gamma, \lambda \geq 0$ may be represented by 
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma \lambda I & -\gamma \lambda I \hl I &-\gamma I & 0 \\
 I & -2\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{F_1(z_k^1) \\ F_2(z_k^2)}
\end{align*}
```


Figure [1](#fig-dr) visualizes  executions of the Douglas-Rachford algorithm to solve the optimization problem $\min f(\beta) + \norm{\beta}_1$ for a quadratic $f$.

:::{figure} img/dr_trace.webp
:alt: Multiple trajectories of the Douglas-Rachford Algorithm
:name: fig_dr

 *Figure 1:* The optimal solution is the black circle. The visualized curves are the outputs $\{z_k^2\}_{k \in \N}$ starting from random initial conditions $x_0$.
:::

## Convergence Properties


A fixed-point of the algorithm is a tuple $(x^*, w^*, z^*)$ satisfying 
```{math}
\begin{align*}
 (F, G): & & \mat{c}{x^* \hl z^*} &= \mat{c|c}{\Acl & \Bcl \hl \Ccl & \Dcl}   \mat{c}{x^* \hl w^*}, & w^* \in  F(z^*).
\end{align*}
``` 

The algorithm is convergent if for every initial condition $x_0$, there exists a fixed point $(x^*(x_0), w^*(x_0), z^*(x_0))$ such that 
1. Optimality: $\sum_{i=1}^s w^{*i}(x_0) = 0$ 
2. Consensus:  $z^{*1}(x_0) = z^{*2}(x_0) = \ldots = z^{*s}(x_0)$
3. Attractivity:  $\lim_{k\rightarrow \infty} \mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2 = 0$.


The algorithm is linearly  convergent in state with rate $\rho \in (0, 1)$ if there exists a constant $\gamma_0> 0$ with
```{math}
\begin{align*}
 \mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2  \leq \gamma_0 \rho^{k} \norm{x_0 - x^*(x_0)}_2 & & \forall k \in \N, \ x_0.
\end{align*}
``` 

