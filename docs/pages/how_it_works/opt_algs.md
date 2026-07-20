# Optimization Algorithms

:::{caution} 
Under Construction
:::

## Optimization/Inclusion Problems

Several problems in optimization, engineering, and economics can be phrased as finding a zero in a sum of $s$ operators
```{math}
\begin{align*}
 0 \in \sum_{i=1}^s F_i(\beta^*).
\end{align*},
```

Zero-inclusion problems include solution concepts such as optimal points and variational inequalities. If each $F_i$ is the subdifferential of a proper function $f_i$, then 
```{math}
\begin{align*}
 \beta^* & \in \argmin_\beta \sum_{i=1}^s f_i(\beta) &  \text{implies} &  \sum_{i=1}^s 0 \in F_i(\beta^*).
\end{align*}
```
The zero inclusion problem is an optimality principle for the composite minimization problem.


An iterate $(w^*, z^*)$ solves the zero-inclusion problem if it satisfies the following properties
1. Consensus: $z^*_1 = z^*_2 = \ldots = z^*_s$,
2. Zeroness: $\textstyle \sum_{i=1}^s w^*_i = 0$,
3. Consistency: $w^*_i \in F_i(z^*_i)$ for all $i \in 1, \ldots, s$.


The operators $\{F_i\}_{i=1}^s$  will be assembled into a single operator $F$, defined as  
```{math}
\begin{align*}
 \text{graph}(F) &= \{(w_i, z_i)_{i=1}^s \mid w_i \in F_i(z_i)\}, & \text{with} &  & w^* &\in F(z^*). 
\end{align*}
```

The solution $\beta^*$ can then be read a valid $z^*$ as $\beta^* = z^*_1$.


## Optimization Algorithms



An optimization algorithm is a procedure that generates a sequence of iterates $(w_k, z_k)_{k \in \N}.$




<!-- Optimization algorithms can be understood as a dynamical is a sequence $(x_k, w_k, z_k)$ -->
<!-- , and include an internal state $x$ -->


 Many common optimization algorithms can be expressed as the interconnection of the operator $F$ with a linear time invariant system $G$ 
```{math}
\begin{align*}
 (F, G): & & \mat{c}{x_{k+1} \hl z_k} &= \mat{c|c}{\Acl & \Bcl \hl \Ccl & \Dcl}   \mat{c}{x_{k} \hl w_k}, & w_k \in  F(z_k).
\end{align*}
``` 

As an example, the  interconnection representation of the Projected Gradient Descent algorithm with parameter $\gamma \geq 0$ is 
```{math}
\begin{align*}
\beta_{k+1} &= \text{proj}_{\mathcal{Z}}(\beta_k - \gamma \nabla f_1(\beta_k)), \\ 
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma I & -\gamma I \hl I &0 & 0 \\
 I & -\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\nabla f_1(z_k^1) \\ \partial f_2(z_k^2)}
\end{align*}
```


