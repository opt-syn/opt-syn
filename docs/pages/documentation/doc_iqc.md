# Integral Quadratic Constraints

Input-output sequences of the operators satisfy a family of relations defined by {doc}Integral Quadratic Constraints (IQC).


---

tocdepth: 1

---


## IQC Analysis


Sequences $(w, z)$ obeying $w_k \in F(z_k)$ are constrained by $F$'s membership in the operator class $\F$. Valid relations among $(w, z)$ may expressed in the framework of Integral Quadratic Constraints (IQC) if $0 \in F(0)$. This IQC methodology is used to certify the Robust Stability portion of the convergence test, the Regulator Equation requirement must be evaluated separately. 


### IQCs and Valid Relations



An IQC with terminal cost is a tuple $(\Psi, M, X)$, where $M, X$ are symmetric matrices, and $\Psi$ is a linear system
```{math}
\begin{align*}
\Psi: & & \mat{c}{\psi_{k+1} \hl r_k} &= \mat{c|cc}{A_\psi& B_\psi^p & B_\psi^q \hl
C_\psi & D_\psi^p & D_\psi^q} \mat{c}{\psi_k \hl p_k \\ q_k}
\end{align*}
```
The system $\Psi$ is referred to as the filter of the IQC.

A sequence $(p, q)$ obeys the IQC $(\Psi, M, X)$ if the relation 
```{math}
\sum_{t=0}^{T-1} r_k^\top M r_k + \psi_{T}^\top X \psi_T \geq 0
```
holds for all time horizons $T \in \N$, whenever  $\Psi$ is driven by the inputs $(p, q)$ starting from the initial condition $\psi_0 = 0$.

### IQCs for Operators

Every operator class obeys a family of valid relations. Each relation is parameterized by a tuple  $(\Psi, M, X, \ell)$, where  $\ell$ is a  signal transformation matrix 
```{math}
\mat{c}{p_k \\ z_k} = \mat{cc}{\ell_{pq} & \ell_{pw} \\ \ell_{zq} & \ell_{qw}} \mat{c}{q_k \\ w_k},
```
and the signals $(p, q)$ derived from $w_k \in F(z_k)$ satisfy the IQC $(\Psi, M, X)$.


As an example, if an $m$-strongly convex function $f_0$ obeys  $0 \in \partial f_0(0)$ and $0 = f_0(0)$, then sequences $(w, z)$ with $w_k \in \partial f_0(z_k)$ for all $k \in \N$ satisfy a relation $(\Psi, M, X, \ell)$ defined by
```{math}
\begin{align*}
\Psi_1: \quad & \mat{c}{\psi_{k+1} \hl r_k} = \mat{c|cc}{0& 0 & I \hl
\lambda_1 I  & \lambda_0 I & 0\\
0 & 0 & I} \mat{c}{\psi_k \hl p_k \\ q_k},  \\
M_1 &:= \mat{c}{0 & I \\ I & 0}, \qquad X_1 = 0, \qquad \ell_1  := \mat{cc}{0 & I \\ I & -m I}, 
\end{align*}
```
for any scalar coefficients $\lambda_0, \lambda_1$ satisfying 
```{math}
\begin{align*}
\lambda_1 & \leq 0, & \lambda_0 + \lambda_1 > 0.
\end{align*}
```


### Exponential Weighting

Linear convergence can be established by proving boundedness of an exponentially weighted system.

Given a rate $\rho > 0$, the  $\rho$-exponential weighting of a sequence $x$ is $\ov x$, defined as $\ov x_k := \rho^{-k} x_k$ for all $k \in \N$. The exponential weighting of an optimization algorithm $(F_0, G)$ with $0 \in F_0(0)$ is 
```{math}
\begin{align*}
 (F_0, G): & & \mat{c}{\ov{x}_{k+1} \hl \ov{z}_k} &= \mat{c|c}{\rhoi \Acl & \rhoi \Bcl \hl \Ccl & \Dcl}   \mat{c}{\ov{x}_{k} \hl \ov{w}_k}, & \ov{w}_k \in  \rho^{-k} F_0(\rho^k \ov{z}_k).
\end{align*}
``` 

If $\rho \in (0, 1)$, then the  existence of a $\gamma_0 > 0$ ensuring boundedness of the exponentially weighted system then proves linear convergence of the original system:

```{math}
\begin{align*}
 \mav{c}{\ov x_k  \\ \ov w_k  \\ \ov z_k }_2 &  \leq \gamma_0 \norm{\ov x_0}_2, & \text{implies} & & \mav{c}{ x_k  \\  w_k  \\  z_k }_2  & \leq \gamma_0 \rho^{k} \norm{x_0}_2.
\end{align*}
``` 

Valid relations for operators can be adjusted to allow for exponential weightings.  The exponentially weighted sequences $(\ov w, \ov z)$ with  $\ov w_k \in \rho^{-k} F(\rho^k \ov z_k)$ for all $k$  satisfies the one-step relation $(\Psi_1, M_1, X_1, \ell_1 )$, but with coefficients 
```{math}
\begin{align*}
\lambda_1 & \leq 0, & \lambda_0 + \rhoi \lambda_1 > 0.
\end{align*}
```

### Outlook for Optimization


IQC descriptions of uncertainties can be used to validate optimization algorithms.


Let $\rho > 0$ be a convergence rate, $\Lambda(\rho)$ be a constraint set, and  $(\Psi(\lambda), M(\lambda), X(\lambda), \ell)$ be a valid relation for any 
1. $F \in \F$ 
2. pairs $(\ov w, \ov z)$ satisfying $\ov w_k \in \rho^{-k} F(\rho^k \ov z_k)$ for all $k$. 
3. $\lambda \in \Lambda(\rho)$. 


The interconnection of linear systems $(\Psi(\lambda), \ell, G)$ with description
```{math}
\begin{align}
\mat{c}{\ov{x}_{k+1} \hl \ov z_k} &= \mat{c|c}{\rhoi \Acl & \rhoi \Bcl \hl \Ccl & \Dcl}\mat{c}{\ov{x}_k \\ \ov w_k},\\
\mat{c}{\ov p_k \\ \ov z_k} &= \mat{cc}{\ell_{pq} & \ell_{pw} \\ \ell_{zq} & \ell_{qw}} \mat{c}{\ov q_k \\ \ov w_k}, \\
\mat{c}{\psi_{k+1} \hl r_k} &= \mat{c|cc}{A_\psi& B_\psi^p & B_\psi^q \hl
C_\psi(\lambda) & D_\psi^p(\lambda) & D_\psi^q(\lambda)} \mat{c}{\psi_k \hl \ov p_k \\ \ov q_k}.
\end{align}
```
is summarized as the system 
```{math}
G(\lambda) := \mat{c}{\psi_{k+1} \\ \ov{x}_{k+1} \hl \ov r_k} &= \mat{c|c}{A(\blam) & B(\blam) \hl C(\blam)& D(\blam)} \mat{c}{\psi_k \\ \ov{x}_k \hl \ov q_k}.
```

A sufficient condition for  $\ov x$ to be bounded for all $F \in \F$ is if there exists a symmetric matrix  $P$ and coefficients $\blam$ such that
```{math}
\mat{cc}{I & 0 \\
A(\blam) & B(\blam) \\
C(\blam) & D(\blam)}^\top \mat{ccc}{-P & 0 & 0\\
0 & P & 0 \\
0 & 0 & M(\blam)} \mat{cc}{I & 0 \\
A(\blam) & B(\blam) \\
C(\blam) & D(\blam)} \prec 0,
```
```{math}
\begin{align*}
P - \mat{cc}{X & 0 \\ 0 & 0 } &\succ 0, & 
\blam & \in \Lambda(\rho).
\end{align*}
```
Boundedness of the $\rho$- weighted $\ov x$ implies linear convergence at rate $\rho$ of the original state $x$ if $\rho \in (0, 1)$. 
Bisection can be used to minimize the rate $\rho$ in this IQC Analysis problem, thus upper-bounding a worst case linear convergence rate. 
Conservatism can be reduced by consider higher-order valid IQCs.
 
## IQCS for Analysis

The IQC structure used in {{osyn}} for analysis of operators is a tuple $(\Psi, M, X, \ell)$, in which $\Psi$ is structured as 

```{math}
\begin{align}
        \Psi = \mat{cc}{\Psi_1 & 0 \\ 0 & \Psi_2}: & & \mat{c}{\psi^1_{k+1} \\ \psi^2_{k+1} \hl p^\psi_{k} \\ q^\psi_{k}} =  \mat{cc|cc}{A_1 & 0 & B_1 & 0 \\
        0 & A_2 & 0 & B_2 \hl
        C_1 & 0 & D_1 & 0 \\
        0 & C_2 & 0 & D_2 
}  \mat{c}{\psi^1_{k+1} \\ \psi^2_{k+1} \hl p_k \\ q_k}.
\end{align}
```

The analysis IQCs may be tall: the dimensions of $(p^\psi, q^\psi)$ may be larger than the dimensions of $(p, q)$ respectively.

## IQCS for Synthesis

IQC synthesis requires filters in which the dimensions of $(p^\psi, q^\psi)$ and $(p, q)$ are equal. The class of filters used are

```{math}
\begin{align}
        \Psi = \mat{cc}{\hat \Psi_1 & \hat \Psi_3 \\ 0 & \hat \Psi_2}: & & \mat{c}{\psi^1_{k+1} \\ \psi^2_{k+1} \hl p^\psi_{k} \\ q^\psi_{k}} =  \mat{cc|cc}{\hat A_1 & 0 &\hat  B_1 & 0 \\
        0 & \hat A_2 & 0 & \hat B_2 \hl
        \hat C_1 & \hat C_3 & \hat D_1 & \hat D_3 \\
        0 & \hat C_2 & 0 & \hat D_2 
}  \mat{c}{\psi^1_{k+1} \\ \psi^2_{k+1} \hl p_k \\ q_k}.
\end{align}
```

In IQC synthesis, it is required that $\Psi_1, \Psi_2, \Psi_3$ are all stable, and that $\Psi_2$ has a stable inverse. A tall IQC from analysis must be factorized into this square form first. 

## IQC Routines

The iqcs routines are `iqc_loop_split` (analysis) and `iqc_loop_factored` (synthesis).

### Analysis Structure

```{eval-rst}
.. mat:autoclass :: iqc.iqc_loop_split   
    :members:
```

### Synthesis Structure



```{eval-rst}
.. mat:autoclass :: iqc.iqc_loop_factored 
    :members:
```

### Data Container

```{eval-rst}
.. mat:autoclass :: iqc.iqc_data_container 
    :members:
```