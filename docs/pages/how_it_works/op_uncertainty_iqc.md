---

tocdepth: 1

---


# Uncertainty Descriptions


Sequences $(w, z)$ obeying $w_k \in F(z_k)$ are constrained by $F$'s membership in the operator class $\F$. Valid relations among $(w, z)$ may expressed in the framework of Integral Quadratic Constraints (IQC) if $0 \in F(0)$. The restriction 


## IQCs and Valid Relations



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
holds for all time horizons $T \in \N$, whenever $\psi_0 = 0$ and $\Psi$ is driven by the inputs $(p, q)$.


Every operator class obeys a family of valid relations. Each relation is parameterized by an IQC $(\Psi, M, X)$ and a signal transformation matrix $\ell$ as
```{math}
\mat{c}{p_k \\ z_k} = \mat{cc}{\ell_{pq} & \ell_{pw} \\ \ell_{zq} & \ell_{qw}} \mat{c}{q_k \\ w_k}.
```


As an example, if an $m$-strongly convex function $f_0$ obeys  $0 \in \partial f_0(0)$ and $0 = f_0(0)$, then sequences $w_k \in \partial f_0(z_k)$ satisfy an IQC defined by
```{math}
\begin{align*}
\Psi: \quad & \mat{c}{\psi_{k+1} \hl r_k} = \mat{c|cc}{0& 0 & I \hl
\lambda_1 I  & \lambda_0 I & 0\\
0 & 0 & I} \mat{c}{\psi_k \hl p_k \\ q_k},  \\
M &:= \mat{c}{0 & I \\ I & 0}, \qquad X = 0, \ell & := \mat{cc}{0 & I \\ I & -m I}, 
\end{align*}
```
for any scalar coefficients $\lambda_0, \lambda_1$ satisfying 
```{math}
\begin{align*}
\lambda_1 & \leq 0, & \lambda_0 + \lambda_1 > 0.
\end{align*}
```

This system $\Psi$ is an order-1 Zames-Falb filter for the operator $F_0$. 

Valid IQCs for operators can be adjusted to allow for exponential weightings of the operators.  The exponentially weighted sequences $(\ov w, \ov z)$ with  $\ov w_k \in \rho^{-k} F(\rho^k \ov z_k)$ satisfies the  IQC $(\Psi, M, X)$ with inputs $(\ov p, \ov q)$, but with coefficients $\lambda_1 \leq 0, \lambda_0 + \rhoi \lambda_1 > 0$.



## Exponential Weighting


The $\rho$-exponential weighting of a sequence $x$ is $\ov x$, defined as $\ov x_k := \rho^{-k} x_k$ for all $k \in \N$. The exponential weighting of an optimization algorithm $(F_0, G)$ with $0 \in F_0(0)$ is 
```{math}
\begin{align*}
 (F_0, G): & & \mat{c}{\ov{x}_{k+1} \hl \ov{z}_k} &= \mat{c|c}{\rhoi \Acl & \rhoi \Bcl \hl \Ccl & \Dcl}   \mat{c}{\ov{x}_{k} \hl \ov{w}_k}, & \ov{w}_k \in  \rho^{-k} F_0(\rho^k \ov{z}_k).
\end{align*}
``` 

The existence of a $\gamma_0 > 0$ ensuring boundedness of the exponentially weighted system then proves linear convergence of the original system:

```{math}
\begin{align*}
 \mav{c}{\ov x_k  \\ \ov w_k  \\ \ov z_k }_2 &  \leq \gamma_0 \norm{\ov x_0}_2, & \text{implies} & & \mav{c}{ x_k  \\  w_k  \\  z_k }_2  & \leq \gamma_0 \rho^{k} \norm{x_0}_2.
\end{align*}
``` 

