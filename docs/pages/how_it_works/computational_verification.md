# Computational Verification

IQCs with terminal cost can be used to prove stability, performance, and invariance properties of dynamical systems. These stability proofs are based on finding solutions to Linear Matrix Inequalites.

## Boundedness/Stability with IQCs


Let $G_R$ =  $(A_R, B_R, C_R, D_R)$ be a linear system, let $\rho> 0$ be a convergence rate, and  $\mathcal{O}$ be an operator.
The exponentially weighted iterates $(\ov p, \ov q)$ related by $\ov p_k \in \rho^{-k} \mathcal{O}(\rho^{-k} q_k)$ for the operator $\mathcal{O}$. This pair $(\ov p, \ov q)$ satisfies the IQC
  $(\Psi(\blam), M(\blam), X(\blam))$ for any  $\blam \in \Lambda(\rho)$.

The $\rho$-weighted interconnection of $(G_R, \mathcal{O})$ is
```{math}
\begin{align}
\mat{c}{\ov{x}_{k+1} \hl \ov p_k} &= \mat{c|c}{\rhoi A_R & \rhoi B_R \hl C_R & D_R}{\ov{x}_k \\ \ov q_k}, & \ov q_k & \in \rho^{-k} \mathcal{O}(\rho^{k} \ov p_k).
\end{align}
```

The system $G_\Psi := \Psi(\blam) \mat{c}{G_R \\ I}$  has  input $q$, output $r$, and states $(\psi, \ov{x})$. The state-space representation of $G_\Psi$ is 
```{math}
\begin{align}
G_\Psi: & &    \mat{c}{\psi_{k+1} \\ \ov{x}_{k+1} \hl r_k} &= \mat{c|c}{A(\blam) & B(\blam) \hl C(\blam)& D(\blam)} \mat{c}{\psi_k \\ \ov{x}_k \hl q_k} \\
& & \mat{c|c}{A & B \hl C& D} & := \mat{cc|c}{
A_\psi & B_\psi^p C_R & B^p_\psi D_R + B^q_\psi \\
0 & \rhoi A_R & \rhoi B_R \hl
C_\psi(\blam) & D^p_\psi(\blam) C_R & D^\psi_p(\blam) D_R + D^\psi_q(\blam)}
\end{align}
```

If there exists a quadratic function $V(\psi, \ov{x})$ and a tolerances $\epsilon \geq 0, \epsilon_q > 0$ such that for all trajectories $(\psi_k, \ov{x}_k)_{k \in \N}$ of $(G_\Psi, \mathcal{O})$ with $\psi_0=0$ we have
```{math}
\begin{align}
V(\psi_{k+1}, \ov{x}_{k+1}) - V(\psi_k, \ov{x}_k) &\leq -r_k^\top M r_k - \epsilon \norm{q_k}^2_2, \\
V(\psi_k, \ov{x}_k) - \psi_k^\top X \psi_k & \geq 0,
\end{align}
```
then we observe for all time horizons $T \in \N$ that 
```{math}
\sum_{k=0}^{T-1} V(\psi_{k+1}, \ov{x}_{k+1}) - V(\psi_{k}, \ov{x}_k) &\leq - \sum_{k=0}^{T-1} r_k^\top  M r_k - \epsilon \norm{q_k}_2^2, \\
V(\psi_{T}, \ov{x}_T) - V(\psi_0, \ov{x}_0) &\leq -\sum_{k=0}^{T-1} r_k^\top  M r_k - \epsilon \norm{q_k}_2^2
```

Adding the valid IQC relation for $\mathcal{O}$ to the above inequality  yields

```{math}
\textstyle \sum_{k=0}^{T-1} V(\psi_{k+1}, \ov{x}_{k+1}) - V(\psi_{k}, \ov{x}_k) &\leq - \sum_{k=0}^{T-1} r_k^\top  M r_k - \epsilon \norm{q_k}_2^2, \\
V(\psi_{T}, \ov{x}_T) - V(\psi_0, \ov{x}_0) &\leq  - \epsilon \textstyle\sum_{k=0}^{T-1} \norm{q_k}_2^2 + \psi_{T}^\top X \psi_T.
```

This relation ensures boundedness of $(\psi, \ov{x})$:
```{math}
0 \leq V(\psi_{T}, \ov{x}_T) - \psi_{T}^\top X \psi_T  \leq  V(\psi_0, \ov{x}_0)- \epsilon \textstyle \sum_{k=0}^{T-1} \norm{q_k}_2^2
```

If there exists $\epsilon_M, \epsilon_X > 0$ such that 
```{math}
\begin{align}
V(\psi_{k+1}, \ov{x}_{k+1}) - V(\psi_k, \ov{x}_k) &\leq -r_k^\top M r_k - \epsilon \norm{q_k}^2_2 - \epsilon_M (\norm{\psi_k}_2^2 + \norm{\ov{x}_k}_2^2), \\ 
V(\psi_k, \ov{x}_k) - \psi_k^\top X \psi_k &\geq \epsilon_X (\norm{\psi_k}_2^2 + \norm{\ov{x}_k}^2_2), 
\end{align}
```
then summing from $0$ to $T-1$ yields the relations

```{math}
\epsilon_X (\norm{\psi_k}_2^2 + \norm{\ov{x}_k}^2_2) &\leq V(\psi_{T}, \ov{x}_T) - \psi_{T}^\top X \psi_T, \\  V(\psi_{T}, \ov{x}_T) - \psi_{T}^\top X \psi_T & \leq  V(\psi_0, \ov{x}_0)-  \textstyle \sum_{k=0}^{T-1} \left(\epsilon \norm{q_k}_2^2 + \epsilon_M \norm{\psi_k}_2^2  + \epsilon_M \norm{\ov{x}_k}_2^2\right),
```

This then proves asymptotic convergence $\lim_{k\rightarrow \infty} (\psi_k, \ov{x}_k) = 0$.


<!-- The boundedness relation can be expanded using the quadratic $V$ as semidefinite constraints
```{math}
0 \preceq \mat{c}{\psi_T \\ \ov{x}_T}^\top \mat{cc}{V_{\psi \ov{x}} - X & V_{\psi \ov{x}} \\ V_{\psi_x}^\top & V_{xx}} \mat{c}{\psi_T \\ \ov{x}_T}  &\leq  V(\psi_0, \ov{x}_0)- \epsilon \sum_{k=0}^{T-1} \norm{q_k}_2^2
``` -->


<!-- If the boundedness inequality holds strictly when starting at any $\ov{x}_0$, then there exists an $\epsilon_x > 0$ -->

## IQC Analysis (General)


The quadratic function $V$ is expanded into 
```{math}
V(\ov{x}, \psi) = \mat{c}{\psi \\ \ov{x}}^\top \mat{cc}{P_{\psi \psi} & P_{\psi \ov{x}} \\ P_{\psi \ov{x}}^\top & P_{xx}} \mat{c}{\psi \\ \ov{x}}.
```

The IQC Analysis to verify asymptotic convergence of the exponentially weighted interconnection is to find a matrix $P$ and coefficients $\blam \in \Lambda(\rho)$ such that
```{math}
\begin{align}
\mat{cc}{I & 0 \\
A(\blam) & B(\blam) \\
C(\blam) & D(\blam)}^\top \mat{ccc}{-P & 0 & 0\\
0 & P & 0 \\
0 & 0 & M(\blam)} \mat{cc}{I & 0 \\
A(\blam) & B(\blam) \\
C(\blam) & D(\blam)} \prec 0
\end{align}
```

```{math}
\begin{align}
\mat{cc}{P_{\psi \psi} - X(\blam) & P_{\psi \ov{x}} \\ P_{\psi \ov{x}}^\top & P_{xx}} & \succ 0
\end{align}
```

Bisection can be used to find a minimal $\rho$ certifying asymptotic convergence.

The analysis program is convex in $\blam$ for fixed $\rho$ if 
- $\Lambda(\rho)$ is convex for each $\rho$,
- either $\Psi(\blam)$ is $\blam$-independent, or $M(\blam)$ is $\blam$-independent. 

Parameterizing $A_\Psi(\blam), B^p_\psi(\blam), B^q_\psi(\blam),$ will cause a bilinearity  between the variables  $\blam$ and $P$. 

## Outlook for Optimization


The operator class $\F$ satisfies a family of IQCs $(\Psi(\blam), M(\blam), X(\blam))$ with respect to signal transformation matrix $\ell$. These $(\Psi(\blam), M(\blam), X(\blam), \ell)$ terms are assembled from the individual operators $\{F_i\}_{i=1}^s$.
 The worst-case linear convergence rate for any algorithm $(F, G)$ with $F \in \F$ can be upper-bounded by solving an IQC analysis problem on the following plant $G_\Psi$  with input $\ov q$ and output $r$: 
```{math}
\begin{align}
\mat{c}{\ov{x}_{k+1} \hl \ov z_k} &= \mat{c|c}{\rhoi \Acl & \rhoi \Bcl \hl \Ccl & \Dcl}\mat{c}{\ov{x}_k \\ \ov w_k},\\
\mat{c}{\ov p_k \\ \ov z_k} &= \mat{cc}{\ell_{pq} & \ell_{pw} \\ \ell_{zq} & \ell_{qw}} \mat{c}{\ov q_k \\ \ov w_k}, \\
\mat{c}{\psi_{k+1} \hl r_k} &= \mat{c|cc}{A_\psi& B_\psi^p & B_\psi^q \hl
C_\psi & D_\psi^p & D_\psi^q} \mat{c}{\psi_k \hl \ov p_k \\ \ov q_k}.
\end{align}
```

Conservatism can be reduced by increasing the filter order $\Psi$ when considering valid relations.
