# Convergence Conditions (Analysis)



The individual operators $F_i$ in the inclusion problem may be contained in known operator classes. Instances of these operator classes include maximally monotone set-valued maps, subgradients of strongly convex functions, gradients of smooth functions, and classes arising from intersections of these properties. These operator classes specify families of optimization problems.


The algorithm $(F,G)$ must achieve the desired convergence properties for all $F$ arising from these operator classes.




Our linear convergence conditions are presented for a set of operators $\F$ satisfying
1. There is a known linear map $\Delta$ with $\Delta \in \F$,
2. Each $F \in \F$ has a unique pair $(\beta^*, w^*)$ solving the inclusion problem.


Uniqueness requires that the map $z \mapsto 0 $ is not in $\F$. Convergence must hold for all $(F, G)$ with $F \in \F.$


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






Sufficient conditions for a well-posed algorithm $(F, G)$ to be linearly convergent at rate $\rho \in (0, 1)$  for all $F \in \F$ are if
1. *Robust Stability:* For all $x_0$ and $F \in \F$ with $0 \in F(0)$, we have $\lim_{k\rightarrow \infty}  \rho^{-2k} \norm{x_k}_2^2 \rightarrow 0$ for the system
```{math}
\begin{align}
\mat{c}{ x_{k+1} \hl {z}_k} &= \mat{c|c}{ \Acl &  \Bcl \hl \Ccl & \Dcl} \mat{c}{ x_{k} \hl {w}_k}, & {w}_k \in  F( {z_k}),
\end{align}
```
2. Solvability of *Regulator Equation:*  There is a solution $\mathbf{X}^*$ to the linear system of equations
```{math}
\begin{align}
\mat{c}{\mathbf{X}^* \\ 0} = \mat{c|cc}{\Acl & 0 & \Bcl N \hl 
\Ccl & \1 \otimes I & \Dcl N} \mat{c}{\mathbf{X}^* \\ I}.
\end{align}
```

## Interpretation

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