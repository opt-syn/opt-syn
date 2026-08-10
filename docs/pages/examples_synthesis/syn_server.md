# Remote Quadratic Programming 

This example considers a strongly convex  optimization problem 
```{math}
\beta^* \in \argmin_{\norm{\beta}_1 \leq \tau} \sum_{i=1}^{5} f_i(\beta),
```
where each function $f_i$ is a quadratic in $S_{m_i, L_i}$.


The overall composite optimization problem is described using the constants
:::{list-table}
:header-rows: 1
* - Oracle
  - $m$
  - $L$
* - 1
  - 0
  - 5
* - 2
  - 1
  - 2
* - 3
  - -2
  - 1
* - 4
  - 1
  - 6
* - 5
  - 1
  - 1
* - 6
  - 0
  - $\infty$
:::

The Synthesized algorithm to solve this problem must have a block-lower-triangular $\Dcl$ matrix. Synthesis of an algorithm with convergence rate $\rho < 0.7679$ is accomplished by using the code
```{literalinclude} ../../../examples/synthesis/syn_server.m
:linenos: true
:caption: Composite Quadratic Programming
:language: matlab
:lines: 1-16
```


Next, the same optimization  problem must be solved in a remote setting. The network dynamics separating the oracles $\partial f$ include time-delays and channel memory. The following code generates an algorithm with $\rho < 0.9341$ in the networked setting.

```{literalinclude} ../../../examples/synthesis/syn_server_channel.m
:linenos: true
:caption: Composite Quadratic Programming with Network Dynamics
:language: matlab
:lines: 1-44
```