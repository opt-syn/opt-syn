# Optimization Algorithms

## Optimization/Inclusion Problems

<!-- Several problems in optimization, engineering, and economics can be phrased as finding a zero in a sum of $s$ operators -->


A zero-inclusion problem tries to find a point $\beta^*$ in the zero in the sum of operators.
```{math}
\begin{align*}
 0 \in \sum_{i=1}^s F_i(\beta^*).
\end{align*}
```

Zero-inclusion problems include solution concepts such as optimization problems and variational inequalities. If each $F_i$ is the subdifferential of a proper function $f_i$, then 
```{math}
\begin{align*}
 \beta^* & \in \argmin_\beta \sum_{i=1}^s f_i(\beta) &  \text{implies} & & 0 & \in \sum_{i=1}^s  F_i(\beta^*).
\end{align*}
```
The zero inclusion problem in $\partial f$ is an optimality principle for the optimization problem.


An iterate $(w^*, z^*)$ solves the zero-inclusion problem if it satisfies 
1. Consensus: $z^*_1 = z^*_2 = \ldots = z^*_s$,
2. Zeroness: $\textstyle \sum_{i=1}^s w^*_i = 0$,
3. Consistency: $w^*_i \in F_i(z^*_i)$ for all $i \in 1, \ldots, s$.

The solution $\beta^*$ can then be read from a valid $z^*$ as $\beta^* = z^*_1$.




## Optimization Algorithms



An optimization algorithm is a procedure that generates a sequence of iterates $(w_k, z_k)_{k \in \N}.$




<!-- Optimization algorithms can be understood as a dynamical is a sequence $(x_k, w_k, z_k)$ -->
<!-- , and include an internal state $x$ -->


 Many common optimization algorithms can be expressed as the interconnection of operators and a linear system. As an example, the gradient descent/forward-step method with stepsize $\gamma > 0$ may be represented by 
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