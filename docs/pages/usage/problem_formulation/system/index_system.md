# System


A System comprises main three attributes:
```matlab
sys = opt_system(Operator_Class, Network, Controller)
```

## Operators 

The `Operator_Class`  argument in the System is an $s$-length  cell array of {doc}`operators classes <../../../documentation/operators/doc_operators>`. 


```{list-table}
:header-rows: 1
* - Operator Class   
  - Class Name
  - Description
* - Set-Valued Maps
  - {class}`op_gen` 
  - monotonicity, cocoercivity, Lipschitzness, Inverse Lipschitzness
* - Subdifferentials
  - {class}`op_sml`
  - Subdifferentials of $S_{m, L}$ with $-\infty < m \leq L \leq \infty$
* - Proper, Closed, Convex
  - {class}`op_pcc`
  - subdifferentials of $S_{0, \infty}$ (e.g. indicator functions)
* - Quadratics
  - {class}`op_quad`
  - Gradients of quadratics  in $S_{m, L}$
```

## Network and Controller


Algorithms to solve inclusion problems $0 \in \sum_{i=1}^s F(\beta^*)$ are modeled using a {doc}`Generalized Plant <../../documentation/plants/doc_genplant>` framework. The signals in this interconnection are:

:::{list-table} 
:header-rows: 1
:widths: 2 8 2 8 2 8 
*   - 
    - State
    -
    - Plant Input
    -
    - Plant Output
*   - $x^N$
    - network
    - $z$
    - input to operators
    - $w$ 
    - output from operators
*   - $x^c$
    - controller
    - $z_p$
    - performance output
    - $w_p$
    - performance input
*   - 
    - 
    - $y$
    - output to controller
    - $u$
    - input from controller
:::


The algorithmic interconnection to solve the inclusion problem $0 \in \sum_{i=1}^s F(\beta^*)$ is modeled as
```{math}
\begin{align*}
\text{Operator}: & & w_k & \in F(z_k), \\
\text{Network}: & & \mat{c}{x^N_{k+1} \hl z_k \\ z_{p k} \\ y_k} &= \mat{c|cc}{A & B_z & B_{z_p} &  B_u \hl 
C_z & D_{zd} & D_{z w_p} &  D_{zu} \\
C_{z_p} & D_{z_p d} & D_{z_p w_p} &  D_{z_p u} \\
C_y & D_{yd} & D_{y w_p} & D_{yu}} \mat{c}{x_k^N \hl w_k \\ w_{pk} \\ u_k}, \\
\text{Controller}: & & \mat{c}{x^c_{k+1} \\ u_k} &= \mat{c|c}{\Ac & \Bc \hl \Cc & \Dc } \mat{c}{x^c_k \\ y_k}.
\end{align*}
```


The network and controller are 


## LASSO Example

The operator class for a composite optimization problem 
```{math}
\beta^* \in \argmin_\beta  f_1(\beta) + \mathbf{I}_{\norm{\cdot}_1 \leq \tau}(\beta),
```
with $f_1 \in S_{1, 10}$ can be specified using
```matlab
op1 = op_sml(1, 10);
op2 = op_pcc();
Operator_Class = {op1, op2};
```

## Extensions




```{toctree}
:maxdepth: 1
:hidden:
Operators <operators>
Dynamical Systems <dynamical_systems/index_systems>
Bind <bind>
Tracking <tracking>
```