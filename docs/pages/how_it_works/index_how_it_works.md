# How it Works

{{osyn}} analyzes and synthesizes optimization algorithms using concepts from robust control. This page summarizes optimization/inclusion problems, convergence conditions, and the mathematical certificates provided by {{osyn}}.
 
 
An optimization problem tries to find a point $\beta^*$ satisfying 
```{math}
\beta^*  \in \argmin_\beta \sum_{i=1}^s f_i(\beta).
```

An inclusion problem is a more general concept than an optimization problem. Inclusion problems try to find a point $\beta^*$ in the zero in the sum of operators.
```{math}
\begin{align*}
 0 \in \sum_{i=1}^s F_i(\beta^*).
\end{align*}
```

A solution (fixed-point) to a zero-inclusion problem is a pair $(\beta^*, w^*)$ with 

```{math}
\begin{align*}
 0 \in \sum_{i=1}^s w^{*,i}, \qquad w^{*,i} \in F_i(\beta^*).
\end{align*}
```

Zero-inclusion problems include solution concepts such as variational inequalities and constrained optimization problems. In particular, if each $F_i$ is the subdifferential of a proper function $f_i$, then 
```{math}
\begin{align*}
 \beta^* & \in \argmin_\beta \sum_{i=1}^s f_i(\beta) &  \text{implies} & & 0 & \in \sum_{i=1}^s  F_i(\beta^*).
\end{align*}
```
The zero inclusion problem in $\partial f$ is a necessary optimality principle for the optimization problem.




## Optimization Algorithms



An optimization algorithm is a procedure that generates a sequence of iterates $(w_k, z_k)_{k \in \N}$ satisfying $w^i_k \in F^i(z^i_k)$.


 Many common optimization algorithms can be expressed as the interconnection of operators and linear systems. As an example, the gradient descent/forward-step method with stepsize $\gamma > 0$ may be represented by 
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k} &= \mat{c|cc}{I & -\gamma \lambda I \hl I & 0 }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{F_1(z_k^1)},
\end{align*}
```
 and the   Douglas-Rachford algorithm {footcite}`douglas1956numerical`   with parameters $\gamma, \lambda \geq 0$ may be represented by 
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma \lambda I & -\gamma \lambda I \hl I &-\gamma I & 0 \\
 I & -2\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{F_1(z_k^1) \\ F_2(z_k^2)}
\end{align*}
```


Figure [1](#fig-dr) visualizes  executions of the Douglas-Rachford algorithm to solve the optimization problem $\min f(\beta) + \norm{\beta}_1$ for a quadratic $f$.

:::{figure} img/dr_trace.webp
:alt: Multiple trajectories of the Douglas-Rachford Algorithm
:name: fig_dr

 *Figure 1:* The optimal solution is the black circle. The visualized curves are the outputs $\{z_k^2\}_{k \in \N}$ starting from random initial conditions $x_0$.
:::

## Convergence Properties


The algorithm is well-posed if the trajectory $(x_k, w_k, z_k)_{k \in \N}$ is unique for all initial conditions $x_0$.


A fixed-point of the algorithm is a tuple $(x^*, w^*, z^*)$ satisfying 
```{math}
\begin{align*}
 (F, G): & & \mat{c}{x^* \hl z^*} &= \mat{c|c}{\Acl & \Bcl \hl \Ccl & \Dcl}   \mat{c}{x^* \hl w^*}, & w^* \in  F(z^*).
\end{align*}
``` 

The algorithm is convergent if for every initial condition $x_0$, there exists a fixed point $(x^*(x_0), w^*(x_0), z^*(x_0))$ such that 
1. Optimality: $\sum_{i=1}^s w^{*,i}(x_0) = 0$ 
2. Consensus:  $z^{*1}(x_0) = z^{*2}(x_0) = \ldots = z^{*s}(x_0)$
3. Attractivity:  $\lim_{k\rightarrow \infty} \mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2 = 0$.


It is linearly convergent  with rate $\rho \in (0, 1)$ if there exists a constant $\gamma_0> 0$ with
```{math}
\begin{align*}
 \mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2  \leq \gamma_0 \rho^{k} \norm{x_0 - x^*(x_0)}_2 & & \forall k \in \N, \ x_0.
\end{align*}
``` 

## Checking Convergence

{{osyn}} certifies linear convergence of well-posed algorithms by checking two conditions:
1. *Robust Stability:* If the algorithm is well-posed and $0 \in F^i(0)$ for all operators $F^i$,  then for any initial state $x_0$, all trajectories of 
```{math}
\begin{align}
\mat{c}{ x_{k+1} \hl {z}_k} &= \mat{c|c}{ \Acl &  \Bcl \hl \Ccl & \Dcl} \mat{c}{ x_{k} \hl {w}_k}, & {w}_k \in  F( {z_k}),
\end{align}
```
satisfy  $\lim_{k \rightarrow \infty} \rho^{-k} x_k = 0$.


2. Solvability of *Regulator Equation:*  For any pair $(\beta^*, w^*)$ with $\sum_{i=1}^s w^{*,i} = 0$, there exists a state $x^*$ satisfying 
```{math}
\begin{align}
\mat{c}{x^* \\ 0} = \mat{c|cc}{\Acl & 0 & \Bcl \hl 
\Ccl & \1 \otimes I & \Dcl} \mat{c}{x^*\\ \beta^*, w^*}.
\end{align}
```

The Robust Stability criterion is an intensive dynamical test, and will be verified using Integral Quadratic Constraints and Linear Matrix Inequality methods. In contrast, the  Regulator Equation can be easily checked by solving a linear system of  equations. Uniqueness of the state $x^*$ is provided by detectability of $(\Acl, \Ccl)$.

The Regulator Equation requirement is independent of the specific operators in $F$. Robust Stability is verified for classes of operators (e.g. $F_2$ is maximal monotone).


Fullfillment of the Regulator Equation and the {{osyn}}-verified Robust Stability requrirements imply that the algorithm is a fixed-point encoding {footcite}`ryu2020uniqueness`: every fixed point of the algorithm is a fixed point of the inclusion problem.

## Networked Setting

Synthesis is posed in terms of a Network separating the oracle $F$ to the controller (you). These networks can model time-delays, cross-talk, channel memory, and other phenomena. The default case of no network dynamics (direct connection to the oracle $F$) is 

```{math}
\mat{c}{
        z_k \\ y_k
    } &= \mat{cc}{0 & I \\ I & 0} \mat{c}{
        w_k \\ u_k
    }.
```

A more general network can be modeled as a linear system. The interconnection between the network and a controller forms the algorithm that interfaces the operator $F$. 

The Regulator Equation condition in the networked setting can be expanded into 

2. *Regulator Equation*:  For any $(\beta^*, w^*)$ with $\sum_{i=1}^s w^{*,i} = 0$, there exists a 
```{math}
\begin{align}    
     \text{Network}: & &    \mat{c|cc:c}{A & 0 & B_w & B_y \hl
        C_z  & \1 & D_{zw} &  D_{zy} \hdl
        C_y & 0 & D_{yz} & D_{yu}} \mat{c}{\Pi\hl I\hdl \Gamma}&=\mat{c}{\Pi \hl 0 \hdl \Phi}, \label{eq:nominal_regulation_control_sys_plant}  \\
\text{Controller}: & & \mat{c|c}{\Ac&\Bc\hl \Cc&\Dc}
\mat{c}{\Theta\\\Phi}  &= \mat{c}{\Theta \\ \Gamma}. \label{eq:nominal_regulation_control_sys_control}    
\end{align}
```

If there does not exist a $(\Pi, \Gamma, \Phi)$ triple satisfying the top equation, then a convergent optimization algorithm cannot be found.

Under the Convergence and Regulator Equation conditions, convergence in all signals is achieved as
```{math}
\begin{align}
    \lim_{k \rightarrow \infty}
    \mat{c}{x_k^N \\ x^c_k \hl  y_k \\u_k} & = \mat{c}{\Pi \\ \Theta \hl \Phi \\ \Gamma} \mat{c}{-\beta^* \\{w}^{*,1} \\ w^{*,2} \\ \vdots \\ w^{*, s-1}}, &   \lim_{k \rightarrow \infty}
    z_k & =  z^* = \mat{c}{\1 \otimes \beta^*}.     
\end{align}
```

Controllers are formed by the interconnection of an internal model {footcite}`francis1976internal` and a designed subcontroller. The internal model is based on the solutions $(\Pi, \Gamma, \Phi)$ of the regulator equations, and ensures that the Regulator Equation requirement is satisfied for *any* subcontroller.

The subcontroller must then ensure that the overall procedure is well-posed and  obeys  Robust Stability condition. Solving for the subcontroller can be accomplished through  IQC synthesis methods {footcite}`veenman2011iqc`.
Synthesis may also involve selecting a solution $(\Pi, \Gamma, \Phi)$ to the regulator equations {footcite}`scherer1997multiobjective`.


## Extensions

The overview is limited to static optimization problems with time-independent memory/stepsize rules. The {doc}`Problem Formulation <../usage/problem_formulation/index_problem_formulation>` section in {doc}`Usage<../usage/index_usage>` documents generalizations to this base construction, including 
- Performance criteria
- Time-varying optimization problems
- Time-varying dynamical systems
- Repeated operator calls
