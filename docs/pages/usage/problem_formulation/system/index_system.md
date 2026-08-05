# Build the System

Algorithms to solve inclusion problems $0 \in \sum_{i=1}^s F(\beta^*)$ are modeled using a {doc}`Generalized Plant <../../../documentation/plants/doc_genplant>` framework. The System (algorithmic interconnection) is specified by the operators $F$, the network, and the controller. 
The System is mathematically described by
```{math}
\begin{align*}
\text{Operator}: & & w_k & \in F(z_k), \\
\text{Network}: & & \mat{c}{x^N_{k+1} \hl z_k \\ z_{p k} \\ y_k} &= \mat{c|cc}{A & B_z & B_{z_p} &  B_u \hl 
C_z & D_{zd} & D_{z w_p} &  D_{zu} \\
C_{z_p} & D_{z_p d} & D_{z_p w_p} &  D_{z_p u} \\
C_y & D_{yd} & D_{y w_p} & D_{yu}} \mat{c}{x_k^N \hl w_k \\ w_{pk} \\ u_k}, \\
\text{Controller}: & & \mat{c}{x^c_{k+1} \hl u_k} &= \mat{c|c}{\Ac & \Bc \hl \Cc & \Dc } \mat{c}{x^c_k \hl y_k},
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

In Simulation, `Operators{i}` is the specific {doc}`operator <../../../documentation/doc_simulation>` $F_i$ used in the inclusion problem. An operator may implement the following methods:
```{list-table}
:header-rows: 1

* - Evaluation
  - Name
  - Operation
* - Forward
  - {meth}`fw`
  - $z \mapsto F_i(z)$,
* - Backward (with parameter $\mathcal{W} \succ 0$)
  - {meth}`bw`
  - $z \mapsto (I - \mathcal{W} F_i)^{-1} (z)$
* - Function 
  - {meth}`f` 
  - $z \mapsto f_i(z)$
```


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
  - Quadratic function $\frac{1}{2} (x - x^*)^\top M (x - x^*)$
  * - Least Squares
  - {class}`op_sim_lsq`
  - Quadratic function $\frac{1}{2} \norm{Ax - b}^2_2$
* - $L_\infty$ (hard) 
  - {class}`op_sim_box`
  - indicator function of  $L_\infty$ ball 
* - $L_1$ (hard) 
  - {class}`op_sim_l1_hard`
  - indicator function of  $L_1$ ball
* - Linear Quadratic Game
  - {class}`op_sim_lq_game`
  - Pseudogradient of game, with agent payoffs  $f_j = \frac{1}{2} x^\top Q_j x + b^\top x_j + e_j$  
```

### Operators Classes

In Analysis and Synthesis, `Operators{i}` is the {doc}`operator class <../../../documentation/operators/doc_operators>` for which operator $F_i$ is a member. 

The two categories of operator classes are general Set-Valued Maps and Subdifferentials.

#### Set-Valued Maps
A general set-valued map is specified by {class}`op_gen`. Fields of {class}`op_gen` define constraints satisfied by all   $w_1 \in F_i (z_1), w_2 \in F_i(z_2)$.
```{list-table}
:header-rows: 1
* - Field name   
  - Parameter
  - Description
* - `monotone`
  - $\mu \in \R$
  - $\langle w_1 - w_2, z_1 - z_2 \rangle \geq \mu \norm{z_1 - z_2}^2_2$
* - `cocoercive`
  - $\beta > 0$
  - $\langle w_1 - w_2, z_1 - z_2 \rangle \geq \beta \norm{w_1 - w_2}^2_2$
* - `lipschitz`
  - $L > 0$
  - $\norm{w_1 - w_2}_2 \leq  L \norm{z_1 - z_2}_2$
* - `inverse_lipschitz`
  - $R > 0$
  - $\norm{z_1 - z_2}_2 \leq  R \norm{w_1 - w_2}_2$
```

Setting the `monotone` field to  $\mu \in \R$ is a description that $F_i - \mu \ \text{Id}$ is maximal monotone. Strong monotonicity is described by  $\mu>0$, and hypo (weak) monotonicity is described by $\mu < 0$.

#### Subdifferentials

The supported subdifferentials are based on properties of proper, closed, convex (p.c.c.) functions. A p.c.c. function $f$ satisfies the properties
```{list-table}
* - Proper
  - $f(x) > -\infty$ everywhere
* - Closed
  - The set $\{x \mid f(x) \leq \gamma\}$ is closed for all $\gamma \in \R$
* - Convex
  - $f( \alpha x + (1-\alpha)y) \leq \alpha f(x) + (1-\alpha) f(y)$ for all $\alpha \in [0, 1]$ and $(x, y)$.
```
The subdifferential  of a p.c.c. function $f$ is the set 
$$ \partial f(x) = \{g \mid f(x) - f(y) \geq \langle g, x-y \rangle, \ \forall (x, y) \}$$ 

Indicator functions of closed,  convex sets are p.c.c.


Given constants $-\infty < m < L \leq \infty$, the set $S_{m, L}$ is the set of functions such that 
1. $f - \frac{m}{2}\norm{\cdot}_2^2$ is p.c.c
2. $\frac{L}{2} \norm{\cdot} - f$ is p.c.c. if $L < \infty$.

$S_{0, \infty}$ is the set of p.c.c. functions. If $m > 0$, then every $f \in S_{m, \infty}$ is strongly convex. If $0 < m < L < \infty$, then every $f \in S_{m, L}$  has $L$-Lipschitz gradients (is $L$-smooth).

The subdifferential of p.c.c. functions is extended to subdifferentials of functions $f \in S_{m, L}$ by 
\begin{align*}
\partial f(x) := \partial \left(f - \frac{m}{2}\norm{x}_2^2\right) + m x.
\end{align*}


Operators arising from subdifferentials are described using the classes
```{list-table}
* - {class}`op_pcc`
  - Subdifferentials of p.c.c. functions
* - {class}`op_sml(m, L)`
  - Subdifferentials of $S_{m, L}$
* - {class}`op_quad(m, L)`
  - Gradients of quadratics in $S_{m, L}$
```

:::{tip}
{class}`op_pcc` is an alias for {class}`op_sml(0, inf)`. 
:::

## Network and Controller

In Simulation and Analysis, the `Controller` is a discrete-time state space system of type `ss`.  
The `Controller` field is ignored in Synthesis, and can therefore be set to `Controller = []`. 


The declaration `Network = []` is used if there are no network dynamics.
If network dynamics are present, then the  `Network` is described by a {class}`genplant` object (see {doc}`genplant documentation <../../../documentation/plants/doc_genplant>` for more details). The attribute {attr}`P` of a {class}`genplant` 
 is a discrete-time state space system of type [ss](https://www.mathworks.com/help/control/ref/ss.html). 
 
 The {class}`genplant` attributes (`nz`, `nzp`, `ny`, `nw`, `nwp`, `nu`) count dimensions of the respective input and output partitions. 

An example {class}`genplant` declaration for a two-operator problem with one performance input and output channel is
```matlab
D = [0, 0, 1, 1, 0;
     0, 0, 1, 0, 1;
     1, 1, 0, 0, 0;
     1, 0, 0, 0, 0;
     0, 1, 0, 0, 0];

n = struct('nz', 2, 'nzp', 1, 'nw', 2, ...
'nw', 2, 'nwp', 1, 'nu', 2);     
p = genplant(ss(D), n);
```

{class}`genplant` can also be called without the `n` argument, letting the dimensions such as  $p.nw$ be set later.



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
