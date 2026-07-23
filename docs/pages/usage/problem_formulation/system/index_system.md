# System


The algorithmic interconnection to solve the inclusion problem $0 \in \sum_{i=1}^s F(\beta^*)$ is
```{math}
\begin{align*}
\text{Operator}: & & w_k & \in F(z_k), \\
\text{Network}: & & \mat{c}{x^N_{k+1} \hl z_k \\ y_k} &= \mat{c|cc}{A & B_z & B_u \hl 
C_z & D_{zd} & D_{zu} \\
C_y & D_{yd} & D_{yu}} \mat{c}{x_k^N \hl w_k \\ u_k}, \\
\text{Controller}: & & \mat{c}{\xi_{k+1} \\ u_k} &= \mat{c|c}{\Ac & \Bc \hl \Cc & \Dc } \mat{c}{\xi_k \\ y_k}.
\end{align*}
```

The System class specifies this interconnection by 

```matlab
sys = opt_system(Operator, Network, Controller)
```

In the case of no network dynamics, the Network field can be left empty as `Network = []`. 


# Operators 

The `Operator`  argument in the System is an $s$-length  `cell array` of {doc}`operators <../../../documentation/operators/doc_operators>`. Each operator $\{F_i\}_{i=1}^s$ is specified by its residing operator class. 

```{list-table}
:header-rows: 1
* - Operator Class   
  - Class Name
  - Description
* - Set-Valued Maps
  - `op_gen` 
  - monotonicity, cocoercivity, Lipschitzness, Inverse Lipschitzness
* - Subdifferentials
  - `op_sml`
  - Subdifferentials of $S_{m, L}$ with $-\infty < m \leq L \leq \infty$
* - Proper, Closed, Convex
  - `op_pcc`
  - subdifferentials of $S_{0, \infty}$ (e.g. indicator functions of closed convex sets)
* - Quadratics
  - `op_quad`
  - Gradients of quadratics  in $S_{m, L}$
```

The operator class for a composite optimization problem 
```{math}
\beta^* \in  f_1(\beta) + \mathbf{I}_{\norm{\cdot}_1 \leq \tau}(\beta),
```
with $f_1 \in S_{m, L}$ can be specified as 

```matlab
op1 = op_sml(m, L);
op2 = op_pcc();
Operator = {op1, op2};
```



# Network and Controller


# Coordinate Lift

# Extensions




```{toctree}
:maxdepth: 1
:hidden:
Operators <operators>
Dynamical Systems <dynamical_systems/index_systems>
Bind <bind>
Tracking <tracking>
```