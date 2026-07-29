# Convergence Conditions (Analysis)

Our computational linear convergence results are based on a subset of zero-inclusion problems with unique solutions.

## Parameterization

The pair $(w^*, z^*)$ may be parameterized using the consensus matrix {footcite}`upadhyaya2025automated`
```{math}
\begin{align}
N &:= \mat{cc}{I & 0 &  \hdots & 0 \\
0 & I & \hdots & 0 \\
0 & 0 & \hdots & I \\
\vdots & \vdots & \ddots & \vdots \\
-I & -I & -I & -I}, & N^\top z^* & = 0,&  \ w^* & \in \text{span}(N).
\end{align}
```

The solution to the zero-inclusion problem can then be described by a pair $(\beta^*, \hat{w}^*)$ with
```{math}
\begin{align}
z^* &= \1 \otimes \beta^*, & w^* &= N\hat{w}^*, &  F(\1 \otimes \beta^*) &= N \hat{w}^*. 
\end{align}
```



## Conditions

Our convergence conditions are based on the  operator class $\F$ satisfying
1. There is a known linear map $\Delta$ with $\Delta \in \F$,
2. Each $F \in \F$ has a unique pair $(\beta^*, \hat{w}^*)$ such that $F(\1 \otimes \beta^*) = N \hat{w}^*$.

Uniqueness requires that the map $z \mapsto 0 $ is not in $\F$. 


Sufficient conditions for a well-posed algorithm $(F, G)$ to be linearly convergent at rate $\rho \in (0, 1)$  for all $F \in \F$ are if
1. *Robust Stability:* For all $\ov x_0$ and $F \in \F$ with $0 \in F(0)$, we have $\lim_{k\rightarrow \infty} \ov x_k \rightarrow 0$ for the system
```{math}
\begin{align}
\mat{c}{\ov x_{k+1} \hl \ov{z}_k} &= \mat{c|c}{\rhoi \Acl & \rhoi \Bcl \hl \Ccl & \Dcl} \mat{c}{\ov x_{k} \hl \ov{w}_k}, & \ov{w}_k \in \rho^{-k} F(\rho^k \ov{z_k}),
\end{align}
```
2. Solvability of *Regulator Equation:*  There is a solution $\mathbf{X}^*$ to the linear system of equations
```{math}
\begin{align}
\mat{c}{\mathbf{X}^* \\ 0} = \mat{c|cc}{\Acl & 0 & \Bcl N \hl 
\Ccl & \1 \otimes I & \Dcl N} \mat{c}{\mathbf{X}^* \\ I}.
\end{align}
```

The Regulation Equation requirement implies the existence of an $x^*$ for any valid $(w^*, z^*)$ satisfying the fixed point relation 
```{math}
\begin{align}
\Acl x^* + \Bcl w^* - x^* &= 0, & z^* = \Ccl x^* + \Dcl w^*.
\end{align}
```

If these two conditions hold, then 
```{math}
\begin{align}
\lim_{k \rightarrow \infty} (w_k, z_k) &= (w^*, z^*), & 
\lim_{k \rightarrow \infty} x_k &= X^* \mat{c}{-\beta^* \\ \hat{w^*}}.
\end{align}
```

These sufficient conditions are also  necessary if all eigenvalues of $\Acl + \Bcl (I - \Dcl \Delta)^{-1} \Dcl \Ccl$ have absolute value below 1. Fullfillment of the Robust Stability and Regulator Equation requrirements imply that the algorithm is a fixed-point encoding {footcite}`ryu2020uniqueness`: every fixed point of the algorithm is a fixed point of the inclusion problem.



The Robust Stability criterion is an intensive dynamical test, and will be verified using Integral Quadratic Constraints and Linear Matrix Inequality methods. In contrast, the  Regulator Equation can be easily checked by solving a linear system of  equations.