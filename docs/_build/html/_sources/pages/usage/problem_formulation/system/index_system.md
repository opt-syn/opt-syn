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


Algorithms to solve inclusion problems $0 \in \sum_{i=1}^s F(\beta^*)$ are modeled using a {doc}`Generalized Plant <../../../documentation/plants/doc_genplant>` framework. The signals in this interconnection are:

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


The algorithmic interconnection to solve the inclusion problem $0 \in \sum_{i=1}^s F(\beta^*)$ is described as
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

The `Controller` field is ignored in Synthesis, and can therefore  be set to `Controller = []`. In Analysis, the `Controller` is a discrete-time state space system of type `ss`.  

The declaration `Network = []` is used if there are no network dynamics.

<!-- If network dynamics are present,  -->

If network dynamics are present, then the  `Network` is described by a {class}`genplant` object. A {class}`genplant` has two attribute:
1. {attr}`P`: State space description $(A, B, C, D)$ of the network
2. {attr}`n`: Dimensions of the partitions $[z, z_p, y]$  $[w, w_p, u]$

The field `P` is a discrete-time state space system of type [ss](https://www.mathworks.com/help/control/ref/ss.html) with sample time $T=1$. The attribute `n` is a struct with integer fields (`nz`, `nzp`, `ny`)  for the dimensions of the output partition and integer fields (`nw`, `nwp`, `nu`) for dimensions of the input partition. 

The {doc}`Templates <../../../documentation/plants/doc_templates>` page documents commands to generates common network structures, such as {class}`bridge_channel_delay` to add  time delays before and after each operator $\{F_i\}_{i =1}^s$. 

## Two-Operator Example

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

Systems for Synthesis with  and without network dynamics are 
```matlab
%add 2-step time delays before and after \partial f1
delay2 = bridge_channel_delay([2, 0], [2, 0]);
sys_delay = opt_system(Operator_Class, delay2, []);

%no network dynamics
sys_no_network = opt_system(Operator_Class, [], []);
```

Systems for Analysis of a Projected Gradient Descent algorithm over the same networks are
```matlab
%the controller describing Projected Gradient Descent
gamma = 2/11;

Ac = 1;
Bc = [-gamma, -gamma];
Cc = [1; 1];
Dc = [0, 0; 
     -gamma, -gamma];
Ts = 1; %sample time     

K = ss(Ac, Bc, Cc, Dc, 1);

sys_pgd_delay = opt_system(Operator_Class, delay2, K);
sys_pgd_no_network = opt_system(Operator_Class, [], K);
```



## Extensions


The System descripition can be extended in three main capacities:


```{toctree}
:maxdepth: 1
Repeated Operator Evaluations <bind>
Time-Varying Optimal Solutions <tracking>
Time-Varying Dynamical Systems <dynamics>
```

These extensions are explored in subsequent sections.
