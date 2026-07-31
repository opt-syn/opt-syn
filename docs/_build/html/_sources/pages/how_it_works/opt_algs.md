# Optimization Algorithms

## Optimization/Inclusion Problems

Several problems in optimization, engineering, and economics can be phrased as finding a zero in a sum of $s$ operators
```{math}
\begin{align*}
 0 \in \sum_{i=1}^s F_i(\beta^*).
\end{align*}
```

Zero-inclusion problems include solution concepts such as optimal points and variational inequalities. If each $F_i$ is the subdifferential of a proper function $f_i$, then 
```{math}
\begin{align*}
 \beta^* & \in \argmin_\beta \sum_{i=1}^s f_i(\beta) &  \text{implies} &  \sum_{i=1}^s 0 & \in F_i(\beta^*).
\end{align*}
```
The zero inclusion problem in $\partial f$ is an optimality principle for the optimization problem in $f$.


An iterate $(w^*, z^*)$ solves the zero-inclusion problem if it satisfies 
1. Consensus: $z^*_1 = z^*_2 = \ldots = z^*_s$,
2. Zeroness: $\textstyle \sum_{i=1}^s w^*_i = 0$,
3. Consistency: $w^*_i \in F_i(z^*_i)$ for all $i \in 1, \ldots, s$.

The solution $\beta^*$ can then be read a valid $z^*$ as $\beta^* = z^*_1$.




## Optimization Algorithms



An optimization algorithm is a procedure that generates a sequence of iterates $(w_k, z_k)_{k \in \N}.$




<!-- Optimization algorithms can be understood as a dynamical is a sequence $(x_k, w_k, z_k)$ -->
<!-- , and include an internal state $x$ -->


 Many common optimization algorithms can be expressed as the interconnection of operators and a linear system. The operators $\{F_i\}_{i=1}^s$  will be assembled into a single operator $F$, defined as  
```{math}
\begin{align*}
 \text{graph}(F) &= \{(w_i, z_i)_{i=1}^s \mid w_i \in F_i(z_i)\}, & \text{with} &  & w^* &\in F(z^*). 
\end{align*}
```
The interconnection of the operator $F$ with a linear time invariant system $G$ is 
```{math}
\begin{align*}
 (F, G): & & \mat{c}{x_{k+1} \hl z_k} &= \mat{c|c}{\Acl & \Bcl \hl \Ccl & \Dcl}   \mat{c}{x_{k} \hl w_k}, & w_k \in  F(z_k).
\end{align*}
``` 

As an example, the  interconnection representation of the Douglas-Rachford algorithm with parameters $\gamma, \lambda \geq 0$ is 
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma \lambda I & -\gamma \lambda I \hl I &-\gamma I & 0 \\
 I & -2\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{F_1(z_k^1) \\ F_2(z_k^2)}
\end{align*}
```


Figure [1](#fig-dr) visualizes  executions of the Douglas-Rachford algorithm to solve the optimization problem $\min f_1(\beta) + \norm{\beta}_1$.

:::{figure} img/dr_trace.webp
:alt: Multiple trajectories of the Douglas-Rachford Algorithm
:name: fig_dr

 *Figure 1:* The optimal solution is the black circle. The visualized curves are the outputs $\{z_k^2\}_{k \in \N}$ starting from random initial conditions $x_0$.
:::
<!-- :scale: 50 % -->


## Operator Classes

The individual operators $F_i$ in the inclusion problem may be contained in known operator classes. The linear system $G$ should ensure that the algorithm $(F, G)$ satisfies the previous properties for all $F$ in a desired class of operators  $\F$. 


The two major classes of supported operators are set-valued maps and subdifferentials of functions.


### Set-Valued Maps
 The set-valued maps may optionally possess the properties of 
- maximal monotonicity,
- $\mu$-strong-monotonicity $\mu > 0$,
- $\mu$-hypo-monotonicity with $\mu > 0$,
- $\beta$-cocoercivity with $\beta>0$,
- $L$-Lipschitzness with $L > 0$,
- $R$-Inverse-Lipschitzness with $R > 0$. 


### Subdifferentials

The most basic subdifferential operator class is the subdifferential of a proper, closed, convex (p.c.c.) function. 
The subdifferential class that we use are characterized by constants $-\infty < m < L < \infty$. The class of functions $\mathcal{S}_{m, L}$ is the set where $f - \frac{m}{2}\norm{\cdot}_2^2$ is p.c.c, and $\frac{L}{2} \norm{\cdot} - f$ is p.c.c. if $L < \infty$. 


{{osyn}} supports
- Subdifferentials of $\mathcal{S}_{m, L}$
- Subdifferentials of quadratics in  $\mathcal{S}_{m, L}$