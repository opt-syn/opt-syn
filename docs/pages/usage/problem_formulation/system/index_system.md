# Build the System

Algorithms to solve inclusion problems $0 \in \sum_{i=1}^s F(\beta^*)$ are modeled using a {doc}`Generalized Plant <../../../documentation/plants/doc_genplant>` framework. The System (algorithmic interconnection) is specified by the operators $F$, the network, and the controller. 



It is mathematically described as
```{math}
\begin{align*}
\text{Operator}: & & w_k & \in F(z_k), \\
\text{Network}: & & \mat{c}{x^N_{k+1} \hl z_k \\ z_{p k} \\ y_k} &= \mat{c|cc}{A & B_z & B_{z_p} &  B_u \hl 
C_z & D_{zd} & D_{z w_p} &  D_{zu} \\
C_{z_p} & D_{z_p d} & D_{z_p w_p} &  D_{z_p u} \\
C_y & D_{yd} & D_{y w_p} & D_{yu}} \mat{c}{x_k^N \hl w_k \\ w_{pk} \\ u_k}, \\
\text{Controller}: & & \mat{c}{x^c_{k+1} \\ u_k} &= \mat{c|c}{\Ac & \Bc \hl \Cc & \Dc } \mat{c}{x^c_k \\ y_k},
\end{align*}
```

where the specific signals are 
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


The System is programatically described by an {class}`opt_system` object:
```matlab
sys = opt_system(Operators, Network, Controller)
```




## Operators 

The `Operators`  argument in the System is an $s$-length cell array `{op1, op2, op3, ...}`.


### Operators for Simulation

In Simulation, `Operators{i}` is the specific {doc}`operator <../../../documentation/doc_simulation>` $F_i$ used in the inclusion problem. An operator 
```{list-table}
:header-rows: 1

* - Evaluation
  - Name
  - Operation
* - Forward
  - {meth}`fw`
  - $z \mapsto F_i(z)$,
* - Backward (with parameter $\Dcl_{ii}$)
  - {meth}`bw`
  - $z \mapsto (I - \Dcl_{ii} F_i)^{-1} (z)$
* - Function 
  - {meth}`f` 
  - $z \mapsto f_i(z)$
```
The {meth}`fw` operation must be defined if $\Dcl_{ii}=0$, and the {meth}`bw` operation must be defined if $\Dcl_{ii} \neq 0$. 


Supported operators for simulation include
```{list-table}
:header-rows: 1
* - Operator    
  - Class Name
  - Description
* - Custom
  - {class}`op_sim` 
  - implemented by [anonymous functions](https://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html)
* - Quadratic
  - {class}`op_sim_quad`
  - Quadratic function $\frac{1}{2} x^\top Q x + b^\top x + e$
* - $L_\infty$ (hard) 
  - {class}`op_sim_box`
  - indicator function of  $L_\infty$ ball 
* - $L_1$ (hard) 
  - {class}`op_sim_box`
  - indicator function of  $L_1$ ball
* - Linear Quadratic Game
  - {class}`op_sim_lq_game`
  - Pseudogradient of game, agent payoffs  $\frac{1}{2} x^\top Q_j x + b^\top x_j + e_j$  
```

### Operators Classes

In Analysis and Synthesis, `Operators{i}` is the {doc}`operator class <../../../documentation/operators/doc_operators>` for which operator $F_i$ is a member.


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

In Simulation and Analysis, the `Controller` is a discrete-time state space system of type `ss`.  
The `Controller` field is ignored in Synthesis, and can therefore be set to `Controller = []`. 


The declaration `Network = []` is used if there are no network dynamics.
If network dynamics are present, then the  `Network` is described by a {class}`genplant` object. The attribute {attr}`P` of a genplant 
 is a discrete-time state space system of type [ss](https://www.mathworks.com/help/control/ref/ss.html) with sample time $T=1$. 
 
 The {class}`genplant` attributes (`nz`, `nzp`, `ny`, `nw`, `nwp`, `nu`) count dimensions of the respective input and output partitions. 

The {doc}`Templates <../../../documentation/plants/doc_templates>` page documents commands to generates common network structures. One such network structure is  {class}`bridge_channel_delay`, which adds time delays before and after each operator $\{F_i\}_{i =1}^s$. 




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
