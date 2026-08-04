# Periodic-Orbit Systems


A Periodic-orbit system has a representation
```{math}
\mat{c}{x_{k+1} \\ z_k} = \mat{c|c}{\Acl_k & \Bcl_k \hl \Ccl_k & \Dcl_k } \mat{c}{x_k \\ w_k}.
```

in which there exists an integer $h$ and $h$-periodic $(M_x, M_w, M_z)$ such that for all $k \in \N$, we have
```{math}
 \mat{c|c}{\Acl_k & \Bcl_k \hl \Ccl_k & \Dcl_k }  = \mat{c c}{M_x & 0 \\ 0 & M_y}^{-k} \mat{c|c}{\Acl_{0} & \Bcl_{0} \hl \Ccl_{0} & \Dcl_{0} } \mat{c c}{M_x & 0 \\ 0 & M_u}^k.
```


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

```{eval-rst}
.. mat:autoclass :: system.periodic_orbit.opt_system_periodic_orbit
    :members:    
```


## Regulator

An open periodic-orbit system with disturbance $d$ and regulated error $e$ is 
```{math}
\begin{align}
d_{k+1} &= S_k d_{k}, \,  \\
\mat{c}{x_{k+1} \hl e_k \\ y_k} &= \mat{c|cc}{A_k & B_{k, d} & B_{k, u} \hl 
C_{k, e} & D_{k, ed} & D_{k, eu} \\
C_{k, y} & D_{k, yd} & D_{k, yu}} \mat{c}{x_k \hl d_k \\ u_k}.
\end{align}
```
The regulator equations for this system are to find matrices $(\Pi, \Gamma, \Phi)$ satisfying
```{math}
\begin{align}
\mat{c}{\Pi M_d S \hl 0 \\ \Phi} &= \mat{c|cc}{M_x A_k & M_x B_{k, d} & M_x B_{k, u} \hl
C_{k, e} & D_{k, ed} & D_{k, eu} \\
C_{k, y} & D_{k, yd} & D_{k, yu}} \mat{c}{\Pi \hl I \\ \Gamma}.
\end{align}
```

If these regulator equations fail, then there does not exist a well-posed and convergent optimization algorithm for this network.

```{eval-rst}
.. mat:autoclass :: system.periodic_orbit.regulator_periodic_orbit
    :members:
```

## LMI Analysis


```{eval-rst}
.. mat:autoclass :: system.periodic_orbit.lmi_analysis_periodic_orbit
    :members:
```

## LMI Synthesis

```{eval-rst}
.. mat:autoclass :: system.periodic_orbit.lmi_synthesis_periodic_orbit
    :members:
```

## LMI Synthesis, Reduced-Order Control

LTI systems allow for reduced-order control synthesis 

```{eval-rst}
.. mat:autoclass :: system.periodic_orbit.lmi_synthesis_periodic_orbit_reduced_order
    :members:
```