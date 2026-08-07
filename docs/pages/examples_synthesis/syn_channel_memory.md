#  Channel Memory

This example continues the {doc}`channel memory simulation <../examples_simulation/sim_channel_memory_symmetric>` example.

An optimization problem
```{math}
\beta^* \in \argmin_{\norm{\beta}_1 \leq 100 } f(\beta)
```


must be solved, where memory effects are present in the communication link to and from evaluation of $\partial f$. 
The network effects with forgetting factor $\alpha > 0$ are
```{math}
\begin{align*}
z^1_k &= u_k^1 - \alpha z^1_{k-1},  & y_k^1 &= w_k^1 - \alpha y_{k-1}^1, \\
z^2_k &= u_k^2, & y_k^2 &= w_k^2, & & \forall k \in \N.
\end{align*}
```

Three rounds of Synthesis/Analysis alternation are performed using `order = {1, 1}` and parameter $\alpha = 0.4$. Upper-bounds on the convergence rate $\rho$ over the course of this alternation are
:::{list-table}
:header-rows: 1
* - Round
  - 1
  - 2
  - 3
* - Synthesis
  - 0.6931    
  - 0.6233    
  - 0.5701
* - Analysis
  - 0.6237
  - 0.5703
  - 0.5348
:::


In contrast, the algorithm from {doc}`channel memory simulation <../examples_simulation/sim_channel_memory_symmetric>` is certified as convergent with $\rho < 0.9448$ under the same `order={1, 1}` Analysis method.


```{literalinclude} ../../../examples/synthesis/syn_channel_symmetric.m
:caption: Code for Repeated Synthesis
:language: matlab
:lines:  1-37
```