# Linear  Time Invariant Systems


A Linear Time Invariant system has a representation
```{math}
\mat{c}{x_{k+1} \\ z_k} = \mat{c|c}{\Acl & \Bcl \hl \Ccl & \Dcl } \mat{c}{x_k \\ w_k}.
```

## System

The algorithmic interconnection is 
```{math}
\begin{align*}
w_k & \in F_k(z_k), \\
\mat{c}{x^N_{k+1} \hl z_k \\ y_k} &= \mat{c|cc}{A & B_z & B_u \hl 
C_z & D_{zd} & D_{zu} \\
C_y & D_{yd} & D_{yu}} \mat{c}{x_k^N \hl w_k \\ u_k}, \\
\mat{c}{\xi_{k+1} \\ u_k} &= \mat{c|c}{\Ac & \Bc \hl \Cc & \Dc } \mat{c}{\xi_k \\ y_k}.
\end{align*}
```

```{eval-rst}
.. mat:autoclass :: system.lti.opt_system
    :members:    
```


## Regulator

An open LTI system with disturbance $d$ and regulated error $e$ is 
```{math}
\begin{align}
d_{k+1} &= S d_k, \\
\mat{c}{x_{k+1} \hl e_k \\ y_k} &= \mat{c|cc}{A & B_d & B_u \hl 
C_e & D_{ed} & D_{eu} \\
C_y & D_{yd} & D_{yu}} \mat{c}{x_k \hl d_k \\ u_k}.
\end{align}
```
The regulator equations for this system are to find $(\Pi, \Gamma, \Phi)$ satisfying
```{math}
\begin{align}
\mat{c}{\Pi S \hl 0 \\ \Phi} &= \mat{c|cc}{A & B_d & B_u \hl
C_e & D_{ed} & D_{eu} \\
C_y & D_{yd} & D_{yu}} \mat{c}{\Pi \hl I \\ \Gamma}.
\end{align}
```

If these regulator equations fail, then there does not exist a well-posed and convergent optimization algorithm for this network.

```{eval-rst}
.. mat:autoclass :: system.lti.regulator_lti
    :members:
```

## LMI Analysis


```{eval-rst}
.. mat:autoclass :: system.lti.lmi_analysis_lti
    :members:
```

## LMI Synthesis

```{eval-rst}
.. mat:autoclass :: system.lti.lmi_synthesis_lti
    :members:
```


## LMI Synthesis, Reduced-Order Control

LTI systems allow for reduced-order control synthesis 

```{eval-rst}
.. mat:autoclass :: system.lti.lmi_synthesis_lti_reduced_order
    :members:
```

<!-- ## System
```{eval-rst}
.. mat:class :: system.lti.opt_system   
    :members:
```

## Regulator
```{eval-rst}
.. mat:class :: system.lti.opt_regulator   
    :members:
```


## Analysis LMIs
```{eval-rst}
.. mat:class :: system.lti.lmi_analysis_lti   
    :members:
```

## Synthesis LMIs
```{eval-rst}
.. mat:class :: system.lti.lmi_synthesis_lti   
    :members:
```


## Reduced-Order Synthesis LMIs
```{eval-rst}
.. mat:class :: system.lti.lmi_synthesis_lti_reduced_order   
    :members:
``` -->