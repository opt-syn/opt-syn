# Properties of Algorithms

Desirable properties of optimization algorithms are convergence, performance, and well-posedness. In Analysis, we aim to verify that an algorithm achieves these properties. In Synthesis, we try to generate an algorithm meeting these properties.

 ## Convergence

A fixed-point of the algorithm $(F, G)$ is a tuple $(x^*, w^*, z^*)$ satisfying 
```{math}
\begin{align*}
 (F, G): & & \mat{c}{x^* \hl z^*} &= \mat{c|c}{\Acl & \Bcl \hl \Ccl & \Dcl}   \mat{c}{x^* \hl w^*}, & w^* \in  F(z^*).
\end{align*}
``` 

The algorithm is convergent if for every initial condition $x_0$, there exists a fixed point $(x^*(x_0), w^*(x_0), z^*(x_0))$ such that 
1. $(w^*(x_0) , z^*(x_0))$ solves the zero-inclusion problem
2. $\lim_{k\rightarrow \infty} \mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2 = 0$.


This convergence is linear with rate $\rho \in (0, 1)$ if there exists a constant $\gamma_0> 0$ with
```{math}
\begin{align*}
 \mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2  \leq \gamma_0 \rho^{k} \norm{x_0 - x^*(x_0)}_2 & & \forall k \in \N, \ x_0.
\end{align*}
``` 

## Performance

The algorithm may be required to operate in noisy and non-ideal environments. Robustness can of the algorithm can be quantified by performance metrics. 


As an example, the addition of  noisy subgradient evaluations $\delta w$ to the interconnection of $(F, G)$ is  
\begin{align*}
  \mat{c}{x_{k+1} \hl z_k} &= \mat{c|c}{\Acl & \Bcl \hl \Ccl & \Dcl}   \mat{c}{x_{k} \hl w_k + \delta w_k}, & w_k \in  F(z_k).
\end{align*}

The algorithm has a finite  $\ell_2$-gain if there exists a $\gamma_x$ such that 
```{math}
\begin{align*}
 \lim_{T \rightarrow \infty} \frac{\sum_{k=0}^T \norm{ x_k - x^*(x_0)}_2^2}{\sum_{k=0}^T \norm{\delta w_k}_2^2}  \leq \gamma_w   & & \text{if  x_0 = x^*(x_0)}.
\end{align*}
``` 

More performance criteria are discussed in {doc}`Performance <../usage/problem_formulation/specs>`.

## Well-Posedness

The interconnection $(F, G)$ may be succinctly expressed using the set-valued map $H := (F^{-1} - \Dcl)^{-1}$ as
```{math}
\begin{align*}
 (F, G): & & x_{k+1} &\in \Acl + \Bcl H(\Ccl x_k).
\end{align*}
``` 
\
The interconnection is well-posed if  $H$ is globally defined and continuous. Well-posedness implies that the sequence $(x_k, w_k, z_k)_{k \in \N}$ exists and is unique for every  initial state $x_0$. A well-posed algorithm therefore satisfies the explicit relation $x_{k+1} = \Acl + \Bcl H(\Ccl x_k)$ for all $k$. 


opt-syn only studies well-posed algorithms. Non well-posed algorithms such as subgradient descent and Frank-Wolfe are currently out of scope.