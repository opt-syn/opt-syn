# Time-Varying Dynamical Systems


A trajectory $(x, w, z)$ of a Linear Time-Varying (LTV) dynamical system with matrix representation $(\Acl_k, \Bcl_k, \Ccl_k, \Dcl_k)_{k \in \N}$ obeys the relation
```{math}
\mat{c}{x_{k+1} \\ z_k} = \mat{c|c}{\Acl_k & \Bcl_k \hl \Ccl_k & \Dcl_k } \mat{c}{x_k \\ w_k}, & \forall k \in \N.
```

Linear Time Invariant (LTI) systems are a special case of LTV systems where $(\Acl, \Bcl, \Ccl, \Dcl)$ are constant in time. 
The {class}`opt_system` class in the {doc}`System <index_system>` page involves networks and controllers that are only  Linear Time Invariant (LTI). 

This page documents time-varying dynamical systems that are supported by {{osyn}}.

In all types of dynamical systems listed on this page,  one of the network or controller can be LTI. If both the network and controller are LTI, then {class}`opt_system` should be used instead of the broader time-varying system class.


## Periodic
A linear time-varying  system is periodic if there exists an integer $h$ such that
```{math}
\begin{align*}
 \mat{c|c}{\Acl_k & \Bcl_k \hl \Ccl_k & \Dcl_k }  = \mat{c|c}{\Acl_{k+h} & \Bcl_{k+h} \hl \Ccl_{k+h} & \Dcl_{k+h} } & & \forall k \in \N.
 \end{align*}
```

Systems with periodic networks and controllers can be specified using the command
```matlab
sys_per = opt_system_periodic(Operator_Class, Network, Controller);
```

An $h$-periodic network is stored as an $h$-length cell of {class}`genplant` objects, or as a {class}`genplant_poly` object. 

If the network is $h$-periodic, then the controller can either be LTI (a single `ss`) or $h$-periodic (an $h$-length cell of `ss` objects).
If the network is LTI, then the controller can be an $h$-length cell of `ss` objects.

A periodic system can be lifted into an LTI system by 
```matlab
sys_lti = sys_per.periodic_lift();
```



## Periodic-Orbit

A periodic-orbit linear system is an $h$-periodic linear system 
in which there exists matrices  $(M_x, M_w, M_z)$ such that 
```{math}
M_x^h &= M_x, \ M_w^h = M_w, M_z^h = M_z, \\
 \mat{c|c}{\Acl_k & \Bcl_k \hl \Ccl_k & \Dcl_k }  &= \mat{c c}{M_x & 0 \\ 0 & M_y}^{-k} \mat{c|c}{\Acl_{0} & \Bcl_{0} \hl \Ccl_{0} & \Dcl_{0} } \mat{c c}{M_x & 0 \\ 0 & M_u}^k, & & \forall k \in \N.
```
Cyclic Coordinate-Descent algorithms are instances of periodic-orbit algorithms satisfying this simplified structure.

The {{osyn}} implementation of periodic-orbit systems is simplified. The matrix $M \in \R^{c \times c}$ must satisfy $M = M^\top, M^h = M$, where $c$ is the coordinate dimension of the Kronecker Structure. All other matrices follow from $M$ by Kronecker structure (e.g. $M_x = I_{n_x/c} \otimes M$).

Systems with periodic-orbit  networks and controllers can be specified using the command
```matlab
sys_orbit = opt_system_periodic_orbit(Operator_Class, Network, Controller, M);
```

Periodic-orbit systems can be enumerated into periodic systems and  lifted into LTI systems:
```matlab
sys_per = sys_orbit.export_periodic();
sys_lti = sys_orbit.periodic_lift();
```

## Switched

A switched linear system with $N_s$ modes is described by a collection of $N_s$ LTI systems (modes), and a directed, unweighted adjacency graph $\mathcal{G}$ with $N_s$ vertices. Each subsystem $j$ of the switched system can be represented as 
```{math}
\mat{c}{x_{k+1} \\ z_k} = \mat{c|c}{\Acl_j & \Bcl_j \hl \Ccl_j & \Dcl_j } \mat{c}{x_k \\ w_k}.
```

The function $\theta: \N \rightarrow \{1, \ldots, N_s\}$ chooses the active subsystem at time $k$. A trajectory $(x, w, z, \theta)_{k \in \N}$ of a switched system  satisfies the relations
```{math}
\begin{align*}
& \mat{c}{x_{k+1} \\ z_k} = \mat{c|c}{\Acl_{\theta(k)} & \Bcl_{\theta(k)} \hl \Ccl_{\theta(k)} & \Dcl_{\theta(k)} } \mat{c}{x_k \\ w_k} & & \forall k \in \N, \\
& (\theta(k), \theta(k+1)) \in \text{Edges}(\mathcal{G}).
\end{align*}
```

Switched systems can model network phenomena such as time-varying delays and communication drops. Periodic and periodic-orbit systems are specific instances of switched systems.

A system with switched system  networks and controllers can be specified using the command
```matlab
sys_per = opt_system_switched(Operator_Class, Network, Controller, G);
```

where `G` is the $\{0, 1\}$ adjacency matrix for the graph $\mathcal{G}$. 