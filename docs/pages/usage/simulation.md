# Simulation




Simulation involves evaluating a trajectory of the system starting from an initial state $x_0$. 


:::{seealso}
- {doc}`Simulation Documentation <../documentation/doc_simulation>` for more details about all objects and routines.
- {doc}`Simulation Examples <../examples_simulation/index_example_simulation>` for demonstrations.

:::

## Execution

Simulation of an inclusion problem $0 \in \sum_{i=1}^s F_i(\beta^*)$ with $\beta^* \in \R^d$ is conducted by the {class}`alg_sim` object. 

The arguments to `alg_sim` are the system (`sys`), and  the dimensionality of $\beta$ in simulation (`d`). The output of an algorithm execution for $T$ time steps is obtained through the {meth}`sim` command,
```matlab
simulator = alg_sim(sys, d);
sim_result = simulator.sim(T);
```


By default, algorithm execution will occur with $x_0=0$, and $w_p = 0$. The  `sampler` field of {class}`alg_sim` allows for random generation and external signals. The attributes of `sampler` are:
:::{list-table}
:header-rows: 1
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
:::


The output of `alg_sim.sim(T)` is an {class}`alg_sim_out` object. The fields of `alg_sim_out` include
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

The {class}`alg_plotter` class accepts a result from simulation, and can then plot any signal stored in  `alg_sim_out`. 
A plot of the signals ($w$, $z$, $u$, $y$) is accomplished by performing
```matlab
plt = alg_plotter(sim_result);
fig = plt.plot({"w", "z", "u", "y"});
```

The figure number can be set by an optional second argument to {meth}`plot`
```matlab
plt = alg_plotter(sim_result);
fig1 = plt.plot({"w", "z", "res_w", "res_z"}, 100); %figure number 100
fig2 = plt.plot({"xi", "u", "y"}, 101); %figure number 101
```


Helper functions of {class}`alg_plotter`  include
:::{list-table} 
 :header-rows: 1
*   - Method
    - Description
    - Plotted Signals
*   - {meth}`plot_6`
    - States, oracles, and convergence
    - (`xn`, `w`, `res_w`, `xc`, `z`, `res_z`)
*   - {meth}`plot_6f`
    - States, oracles, convergence, function values
    - (`x`, `w`, `res_w`, `f`, `z`, `res_z`)
*   - {meth}`plot_4`
    -  oracles, convergence
    - (`w`, `res_w`, `z`, `res_z`)
:::


:::{seealso}
{doc}`Plotting <../documentation/doc_plotting>` contains a  full list of helper functions and more details.
:::

## Details of  Execution

The System `sys` executed  by first interconnecting the Network and Controller, and then interconnecting the possibly nonlinear operator $F$. This is mathematically described by 
```{math}
\begin{align*}
\text{Operator}: & & w_k & \in F(z_k), \\
\text{Algorithm}: & & \mat{c}{x_{k+1} \hl z_k \\ z_{p k}} &= \mat{c|cc}{\Acl & \Bcl_z & \Bcl_{z_p} \hl 
\Ccl_z & \Dcl_{zw} & \Dcl_{z w_p}  \\
\Ccl_{z_p} & \Dcl_{z_p w} & \Dcl_{z_p w_p}} \mat{c}{x_k \hl w_k \\ w_{p k}},
\end{align*}
```
Well-posedness requires that the map $H: = (F^{-1} - \Dcl{zw})^{-1}$ is globally defined and continuous.

The closed loop state $x$ is the concatenation $x = [x^N, x^c]$. 


{class}`alg_sim` execution using {meth}`sim`  requires   well-posedness and a block-triangular information structure (closed loop $\Dcl$ matrix). 

Under these conditions, the System can be partitioned as 

\begin{align*}
\text{Operator}: & & w_k & \in F(z_k), \\
\text{Algorithm}: & & \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2 \\ \vdots \\ z_k^s \\ z_{p k}} &= \mat{c|cccc:c}{\Acl & \Bcl_{z,1} & \Bcl_{z,2} & \cdots & \Bcl_{z,s} & \Bcl_{z_p} \hl 
\Ccl_{z 1} & \Dcl_{zw,11} & 0 & \cdots & 0 & \Dcl_{z w_p, 1}  \\
\Ccl_{z 1} & \Dcl_{zw,21} & \Dcl_{zw,22} & \cdots & 0 & \Dcl_{z w_p, 2}  \\
\vdots &  \vdots & \vdots &  \ddots & \vdots & \vdots  \\
\Ccl_{z 1} & \Dcl_{zw, s1} & \Dcl_{zws2} & \cdots & \Dcl_{zw,ss} & \Dcl_{z w_p s}  \\
\Ccl_{z_p} & \Dcl_{z_p w,1} & \Dcl_{z_p w, 2} & \cdots & \Dcl_{z_p w, s} & \Dcl_{z_p w_p}} \mat{c}{x_k \hl w_k^1 \\ w_k^2 \\ \vdots \\ w_k^s \\ w_{p k}},
\end{align*}


<!-- Algorithm simulation assumes that the System forms a well-posed algorithm with a block-lower triangular  -->

The iterative loop for algorithm simulation is to evaluate the equations for each $k \in \N$
```{math}
\begin{align}
w^i_k &= (F_i^{-1} - \Dcl_{zw, ii})^{-1} (\Ccl_i x_k + \textstyle \sum_{j=1}^{i-1} \Dcl_{zw, ij} w^j_k),  & & \forall i \in 1, \ldots, s,\\ 
x_{k+1} &= \Acl x_k + \textstyle \sum_{i=1}^s \Bcl_{wi} w^i_k + \Bcl_{w_p} w_{pk}. \\ 
z_{p, k+1} &= \Ccl_{z_p} x_k + \textstyle \sum_{i=1}^s \Dcl_{z_p w i} w^i_k + \Dcl_{z_p w_p} w_{pk}.
\end{align}
```



The operator $H_i: = (F_i^{-1} - \Dcl_{zw, ii})^{-1}$ can be evaluated using the {class}`op_sim` methods 
:::{list-table}
:header-rows: 1
* - Evaluation
  - Used Method
  - Condition
  - Operation $z_i \mapsto H_i z_i$
* - Explicit
  - {meth}`fw`
  - $\Dcl_{zw, ii} = 0$
  - $z \mapsto F_i(z_i)$,
* - Implicit
  - {meth}`bw`
  - $\Dcl_{zw,ii}$ is invertible
  - $z_i \mapsto \Dcl_{zw, ii}^{-1} [z - (I - \Dcl_{zw, ii} F_i)^{-1}(z_i)]$
:::

<!-- If the backward-evaluation  $(I - \Dcl_{ii} F_i)^{-1}$ is available (such as from a  resolvent/proximal operator), then this algorithm execution is tractable. -->


