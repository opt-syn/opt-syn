# Periodic Systems


A Periodic system has a representation
```{math}
\mat{c}{x_{k+1} \\ z_k} = \mat{c|c}{\Acl_k & \Bcl_k \hl \Ccl_k & \Dcl_k } \mat{c}{x_k \\ w_k}.
```

in which there exists a period $h$ such that
```{math}
 \mat{c|c}{\Acl_k & \Bcl_k \hl \Ccl_k & \Dcl_k }  = \mat{c|c}{\Acl_{k+h} & \Bcl_{k+h} \hl \Ccl_{k+h} & \Dcl_{k+h} } 
```

A periodic system is a switched system restricted to switching in a ring graph.

## System

The algorithmic interconnection for a periodic system is 
```{math}
\begin{align*}
w_k & \in F_k(z_k), \,  \\
 \mat{c}{x^N_{k+1} \hl z_k \\ y_k} &= \mat{c|cc}{A_{k} & B_{k, \,  z} & B_{k, \,  u} \hl 
C_{k, \,  z} & D_{k, \,  zd} & D_{k, \,  zu} \\ C_{k, \,  y} & D_{k, \,  yd} & D_{k, \,  yu}} \mat{c}{x_k^N \hl w_k \\ u_k}, \\
 \mat{c}{\xi_{k+1} \\ y_k} &=  \mat{c|c}{A_{K, k} & B_{K, k} \hl C_{K, k} & D_{K, k} } \mat{c}{\xi_k \\ y_k}
\end{align*}
``` 

 <!-- & -->
<!-- \mat{c}{\xi_{k+1} \\ u_k} &= \mat{c|c}{\Ac_k & \Bc_k \hl \Cc & \Dc_k } \mat{c}{\xi_k \\ y_k} -->

<!-- 
C_{k, \,  y} & D_{k, \,  yd} & D_{k, \,  yu}} \mat{c}{x_k^N \hl z_k \\ u_k}, -->

```{eval-rst}
.. mat:autoclass :: system.periodic.opt_system_periodic
    :members:    
```


## Regulator

An open periodic system with disturbance $d$ and regulated error $e$ is 
```{math}
\begin{align}
d_{k+1} &= S_k d_{k}, \,  \\
\mat{c}{x_{k+1} \hl e_k \\ y_k} &= \mat{c|cc}{A_k & B_{k, d} & B_{k, u} \hl 
C_{k, e} & D_{k, ed} & D_{k, eu} \\
C_{k, y} & D_{k, yd} & D_{k, yu}} \mat{c}{x_k \hl d_k \\ u_k}.
\end{align}
```
The one-step regulator equations for this system are to find $h$-periodic matrices $(\Pi_h, \Gamma_h, \Phi_h)_{h\in \N}$ satisfying
```{math}
\begin{align}
\mat{c}{\Pi_{k+1} S_k \hl 0 \\ \Phi} &= \mat{c|cc}{A_k & B_{k, d} & B_{k, u} \hl
C_{k, e} & D_{k, ed} & D_{k, eu} \\
C_{k, y} & D_{k, yd} & D_{k, yu}} \mat{c}{\Pi_k \hl I \\ \Gamma_k}.
\end{align}
```

If these regulator equations fail, then there does not exist a well-posed and convergent optimization algorithm for this network.

```{eval-rst}
.. mat:autoclass :: system.periodic.regulator_periodic
    :members:
```

## LMI Analysis


```{eval-rst}
.. mat:autoclass :: system.periodic.lmi_analysis_periodic
    :members:
```

## LMI Synthesis

```{eval-rst}
.. mat:autoclass :: system.periodic.lmi_synthesis_periodic
    :members:
```