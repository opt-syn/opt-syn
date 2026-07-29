# Simulation




Simulation involves evaluating a trajectory of the system starting from an initial state $x_0$. 



 See {doc}`Simulation <../documentation/doc_simulation>` for more details about all objects and routines.

## Execution

Simulation an inclusion problem with $\beta \in \R^d$ is conducted by the {class}`alg_sim` class. The arguments to `alg_sim` are the system `sys`, the dimensionality of $\beta$ `d`. The output of an algorithm execution for $T$ time steps is obtained through the {meth}`sim` command,
```matlab
simulator = alg_sim(sys, d);
sim_result = simulator.sim(T);
```


By default, algorithm execution will occur with zero initial condition and a performance input ($x_0=0$, $w_p = 0$). The  `sampler` field of {class}`alg_sim` has allows for random generation and external signals. 
```{list-table}
:header-rows: 1
:caption:
* - Field
  - Description
* - `x0`
  - Initial Condition
* - `wp`
  - Performance Input  
* - `param0`
  - Initial value of problem-dependent parameters
* - `param`  
  - Subsequent value of problem-dependent parameters
```

The output of `alg_sim.sim(T)` is a {class}`alg_sim_out` object. The fields of `alg_sim_out` include
:::{list-table} 
:widths: 2 8 2 8 2 8 
*   - `xn`
    - state of network
    - `z`
    - input to operators
    - `w` 
    - output from operators
*   - `xi`
    - state of controller
    - `zp`
    - performance output
    - `wp`
    - performance input
*   - `k`
    - time index
    - `y`
    - output to controller
    - `u`
    - input from controller
*   - `f`
    - function value
    - `res_w`
    - optimality error $\norm{\sum_{i=1}^s w^i_k}_2$
    - `res_z`
    - consensus error $\norm{z^i_k - z^i_{\text{average}, k}}_2$
  * - `mode`
    - subsystem for switched systems
    - `param`
    - problem-dependent parameters
    - 
    - 
:::


 Each field of `alg_sim_out` is a numeric or cell array with last dimension indexed by $k \in \{1, \ldots, T\}$.


## Plotting

The {class}`alg_plotter` class accepts a result from simulation.
A plot of the signals "`w`, `z`, `res_w`, `res_z`" is accomplished by performing
```matlab
plt = alg_plotter(sim_result);
fig = plt.plot({"w", "z", "res_w", "res_z"});
```

The figure number can be set by an optional second argument to {meth}`plot`
```matlab
plt = alg_plotter(sim_result);
fig1 = plt.plot({"w", "z", "res_w", "res_z"}, 100); %figure number 100
fig2 = plt.plot({"xn", "u", "y"}, 101); %figure number 101
```

Plottable signals derived from the fields of {class}`alg_sim_out` include `x` (closed-loop state) and `delay` (`mode`-1, used for switched systems).


Given a fixed point 

## Details of  Execution

The System `sys` executed  by first interconnecting the Network and Controller (linear systems), and then interconnecting the operator $F$ (nonlinearity). This is mathematically described by 
```{math}
\begin{align*}
\text{Operator}: & & w_k & \in F(z_k), \\
\text{Algorithm}: & & \mat{c}{x_{k+1} \hl z_k \\ z_{p k}} &= \mat{c|cc}{\Acl & \Bcl_z & \Bcl_{z_p} \hl 
\Ccl_z & \Dcl_{zw} & \Dcl_{z w_p}  \\
\Ccl_{z_p} & \Dcl_{z_p w} & \Dcl_{z_p w_p}} \mat{c}{x_k \hl w_k \\ w_{p k}},
\end{align*}
```
The closed loop state is the concatenation $x = [x^N, x^c]$. 

The System is well-posed if $(F^{-1} - \Dcl)^{-1}$ is invertible, and has a lower-triangular information structure if $\Dcl_{zw}$ is block-lower-triangular (see {doc}`Algorithm Properties <../how_it_works/alg_properties>`). Under these conditions, the System can be partitioned as 

\begin{align*}
\text{Operator}: & & w_k & \in F(z_k), \\
\text{Algorithm}: & & \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2 \\ \vdots \\ z_k^s \\ z_{p k}} &= \mat{c|cccc:c}{\Acl & \Bcl_{z1} & \Bcl_{z2} & \cdots & \Bcl_{zs} & \Bcl_{z_p} \hl 
\Ccl_{z 1} & \Dcl_{zw11} & 0 & \cdots & 0 & \Dcl_{z w_p 1}  \\
\Ccl_{z 1} & \Dcl_{zw21} & \Dcl_{zw22} & \cdots & 0 & \Dcl_{z w_p 2}  \\
\vdots &  \vdots & \vdots &  \ddots & \vdots & \vdots  \\
\Ccl_{z 1} & \Dcl_{zws1} & \Dcl_{zws2} & \cdots & \Dcl_{zwss} & \Dcl_{z w_p s}  \\
\Ccl_{z_p} & \Dcl_{z_p w 1} & \Dcl_{z_p w 2} & \cdots & \Dcl_{z_p w s} & \Dcl_{z_p w_p}} \mat{c}{x_k \hl w_k^1 \\ w_k^2 \\ \vdots \\ w_k^s \\ w_{p k}},
\end{align*}


<!-- Algorithm simulation assumes that the System forms a well-posed algorithm with a block-lower triangular  -->

Algorithm simulation then causally proceeds for each $k \in \N$ as
```{math}
\begin{align}
w^i_k &= (F_i^{-1} - \Dcl_{zw, ii})^{-1} (\Ccl_i x_k + \textstyle \sum_{j=1}^{i-1} \Dcl_{zw ij} w^j_k),  & & \forall i \in 1, \ldots, s,\\ 
x_{k+1} &= \Acl x_k + \textstyle \sum_{i=1}^s \Bcl_{wi} w^i_k + \Bcl_{w_p} w_{pk}. \\ 
z_{p, k+1} &= \Ccl_{z_p} x_k + \textstyle \sum_{i=1}^s \Dcl_{z_p w i} w^i_k + \Dcl_{z_p w_p} w_{pk}.
\end{align}
```

The operator $H_i: = (F_i^{-1} - \Dcl_{zw, ii})^{-1}$ can be evaluated using the {class}`op_sim` methods 
```{list-table}
:header-rows: 1
* - Evaluation
  - Method
  - Condition
  - Operation $H_i$
* - Explicit
  - {meth}`fw`
  - $\Dcl_{zw, ii} = 0$
  - $z \mapsto F_i(z)$,
* - Implicit
  - {meth}`bw`
  - $\Dcl_{ii}$ invertible
  - $z \mapsto \Dcl_{zw ii}^{-1} z - \Dcl_{zw ii}^{-1} (I - \Dcl_{zw ii} F_i)^{-1}(z)$
```

<!-- If the backward-evaluation  $(I - \Dcl_{ii} F_i)^{-1}$ is available (such as from a  resolvent/proximal operator), then this algorithm execution is tractable. -->


