# Synthesis and Networks

Algorithm synthesis aims to design a convergent optimization algorithm obeying the desired {doc}`properties <alg_properties>`. 



## Networked Algorithms

Synthesis is posed in terms of a Network separating the oracle $F$ to the controller (you). These networks can model time-delays, cross-talk, channel memory, and other phenomena.


The default case of no network dynamics (direct connection to the oracle $F$) is 

```{math}
\mat{c}{
        z_k \\ y_k
    } &= \mat{cc}{0 & I \\ I & 0} \mat{c}{
        w_k \\ u_k
    }.
```

The general network $P$ and controller $K$ have the descrptions
```{math}
\begin{align}
    P: & & \mat{c}{
    x^N_{k+1}  \hl z_k \\ y_k
    } &= \mat{c|cc}{ A & B_w & B_u \hl
    C_z & D_{zw} & D_{zu} \\
    C_y & D_{yw} & D_{yu} }\mat{c}{
        x^N_k \hl w_k \\ u_k
    } \\ \\ 
    \Kc : & & \mat{c}{
        \xi_{k+1}  \hl u_k}&=\mat{c|c}{\Ac & \Bc \hl
    \Cc & \Dc}\mat{c}{
        \xi_k \hl y_k
    }.
    \end{align}
```

The interconnection of $(P, K)$ is well-posed if $I - D_{zw} \Dc$ is invertible.


## Convergence Conditions

A subsidiary goal in algorithm synthesis is to ensure that the interconnection  $(F, (P, K))$ satisfies the {doc}`convergence conditions <convergence_conditions>` for optimization algorithms.

The Regulator Equation condition in the networked setting can be expanded into 


2. *Regulator Equation*:  There is a unique solution $(\Pi, \Gamma, \Phi, \Theta)$ to the linear system of equations
```{math}
\begin{align}    
        \mat{c|cc:c}{A & 0 & B_w N & B_y \hl
        C_z  & \1 & D_{zw}N &  D_{zy} \hdl
        C_y & 0 & D_{yz}N & D_{yu}} \mat{c}{\Pi\hl I\hdl \Gamma}&=\mat{c}{\Pi \hl 0 \hdl \Phi}, \label{eq:nominal_regulation_control_sys_plant}  \\
\mat{c|c}{\Ac&\Bc\hl \Cc&\Dc}
\mat{c}{\Theta\\\Phi}  &= \mat{c}{\Theta \\ \Gamma}. \label{eq:nominal_regulation_control_sys_control}    
\end{align}
```

If there does not exist a $(\Pi, \Gamma, \Phi)$ triple satisfying the top equation, then a convergent optimization algorithm cannot be found.

Under the Convergence and Regulator Equation conditions, convergence in all signals is achieved as
```{math}
\begin{align}
    \lim_{k \rightarrow \infty}
    \mat{c}{x_k^N \\ \xi_k \hl  y_k \\u_k} & = \mat{c}{\Pi \\ \Theta \hl \Phi \\ \Gamma} \mat{c}{-\beta^* \\ \hat{w}^*}, &   \lim_{k \rightarrow \infty}
    \mat{c}{z_k \\ w_k } & =  \mat{c}{z^* \\ w^*} = \mat{c}{\1 \otimes \beta^* \\ N \hat{w}^*}.     
\end{align}
```

The fixed-point matrix for the closed-loop system is
```{math}
\begin{align}
\mathbf{X}^* & := \mat{c}{\Pi \\ \Theta}, & \lim_{k \rightarrow \infty} x_k & := \mathbf{X}^* \mat{c}{-\beta^* \\ \hat{w}^*}
\end{align}
```

## Synthesis

The Regulator Equation requirement enforces structure on possible controllers. We choose controllers $\Kc$ as the interconnection of an internal model {footcite}`francis1976internal` and a subcontroller $K_f$. This interconnection with $K = (\text{Model}, K_f)$ involves
```{math}
\begin{align*}
    \text{Model}: \mat{c}{
    v_{k+1} \hl  \tilde{y}_k
    } &= \mat{c|c:cc}{ I & 0 & I & 0 \hl
    -\Gamma & 0 & 0 &  I \\
    \Phi & I & 0 & 0  }\mat{c}{
        v_k \hl y_k \hdl \tilde{u}_k^1 \\ \tilde{u}_k^2
    }, \\
    K_f: \mat{c}{
        \xi^f_{k+1}  \hl \tilde{u}^1_k \\ \tilde{u}^1_k} &=\mat{c|c}{A_f & B_f \hl
    C_{f1} & D_{f1} \\
    C_{f2} & D_{f2}}\mat{c}{
        \xi^f_k \hl \tilde{y}_k 
    }.
    \end{align*}
```

<!-- With this internal model,  -->

Any subcontroller $K_f$ will ensure that the Regulator Equation requirement in $K$ is  satisfied. The subcontroller $K_f$ must then ensure that the interconnection  $(F, (P, \text{Model}, K_f))$ is well-posed and  obeys the convergence and performance requirements. Solving for $K_f$ can be accomplished through  IQC synthesis methods {footcite}`veenman2011iqc`.


Synthesis may also involve selecting a solution $(\Pi, \Gamma, \Phi)$ to the regulator equations {footcite}`scherer1997multiobjective`.

```{footbibliography}
```