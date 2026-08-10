# Switched


A switched linear system with $N_s$ modes is described by a collection of $N_s$ LTI systems (modes), and a directed, unweighted adjacency graph $\mathcal{G}$ with $N_s$ vertices. Each subsystem $j$ of the switched system can be represented as 
```{math}
\mat{c}{x_{k+1} \\ z_k} = \mat{c|c}{\Acl_j & \Bcl_j \hl \Ccl_j & \Dcl_j } \mat{c}{x_k \\ w_k}.
```

The function $\theta: \N \rightarrow 1, \ldots, N_s$ chooses the active subsystem at time $k$. A trajectory $(x, w, z, \theta)_{k \in \N}$ of a switched system  satisfies the relations
```{math}
\begin{align*}
& \mat{c}{x_{k+1} \\ z_k} = \mat{c|c}{\Acl_{\theta(k)} & \Bcl_{\theta(k)} \hl \Ccl_{\theta(k)} & \Dcl_{\theta(k)} } \mat{c}{x_k \\ w_k} & &  \\
& (\theta(k), \theta(k+1)) \in \text{Edges}(\Gs)  & & \forall k \in \N.
\begin{align*}
```


## System

The algorithmic interconnection is 
```{math}
\begin{align*}
w_k & \in F_k(z_k), \,  \\
 \mat{c}{x^N_{k+1} \hl z_k \\ y_k} &= \mat{c|cc}{A_{\theta(k)} & B_{\theta(k), \,  z} & B_{\theta(k), \,  u} \hl 
C_{\theta(k), \,  z} & D_{\theta(k), \,  zd} & D_{\theta(k), \,  zu} \\ C_{\theta(k), \,  y} & D_{\theta(k), \,  yd} & D_{\theta(k), \,  yu}} \mat{c}{x_k^N \hl w_k \\ u_k}, \\
 \mat{c}{\xi_{k+1} \\ y_k} &=  \mat{c|c}{A_{K, \theta(k)} & B_{K, \theta(k)} \hl C_{K, \theta(k)} & D_{K, \theta(k)} } \mat{c}{\xi_k \\ y_k}.
\end{align*}
``` 

 <!-- & -->
<!-- \mat{c}{\xi_{k+1} \\ u_k} &= \mat{c|c}{\Ac_k & \Bc_k \hl \Cc & \Dc_k } \mat{c}{\xi_k \\ y_k} -->

<!-- 
C_{\theta(k), \,  y} & D_{\theta(k), \,  yd} & D_{\theta(k), \,  yu}} \mat{c}{x_k^N \hl z_k \\ u_k}, -->

```{eval-rst}
.. mat:autoclass :: system.switched.opt_system_switched
    :members:    
```


## Regulator

The subsystems for an open switched system with disturbance $d$ and regulated error $e$ may be described as 
```{math}
\begin{align}
d_{k+1} &= S_j d_k, \,  \\
\mat{c}{x_{k+1} \hl e_k \\ y_k} &= \mat{c|cc}{A_j & B_{j, d} & B_{j, u} \hl 
C_{j, e} & D_{j, ed} & D_{j, eu} \\
C_{j, y} & D_{j, yd} & D_{j, yu}} \mat{c}{x_k \hl d_k \\ u_k}.
\end{align}
```
The one-step regulator equations for this system are to find $(\Pi_j, \Gamma_j, \Phi_j)_{j=1}^{N_s}$ satisfying
```{math}
\begin{align}
\mat{c}{\Pi_{j} S_i \hl 0 \\ \Phi_i} &= \mat{c|cc}{A_i & B_{i, d} & B_{i, u} \hl
C_{i, e} & D_{i, ed} & D_{i, eu} \\
C_{i, y} & D_{i, yd} & D_{i, yu}} \mat{c}{\Pi_i \hl I \\ \Gamma_i}.
\end{align}
```
over all arcs $(i, j)$ in the graph $\mathcal{G}$.

Certification of robust stability and feasibility of these regulator equations are sufficient but not necessary to prove convergence. 
If these regulator equations fail, then there may exist a well-posed and convergent optimization algorithm for this network, but {{osyn}} will not find it. 

```{eval-rst}
.. mat:autoclass :: system.switched.regulator_switched
    :members:
```

## LMI Analysis


```{eval-rst}
.. mat:autoclass :: system.switched.lmi_analysis_switched
    :members:
```

## LMI Synthesis

```{eval-rst}
.. mat:autoclass :: system.switched.lmi_synthesis_switched
    :members:
```


