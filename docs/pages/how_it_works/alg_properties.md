# Properties of Algorithms

Desirable properties of optimization algorithms are convergence, performance,  well-posedness, and information structures. In Analysis, we aim to verify that an algorithm achieves these properties. In Synthesis, we try to generate an algorithm meeting these properties. 

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


As an example,  noise $\delta w$ can be introduced into the (sub)-gradient evaluations. The noisy algorithmic loop is 
\begin{align*}
  \mat{c}{x_{k+1} \hl z_k} &= \mat{c|c}{\Acl & \Bcl \hl \Ccl & \Dcl}   \mat{c}{x_{k} \hl w_k + \delta w_k}, & w_k \in  F(z_k).
\end{align*}


The system is input-to-state stable if there exists gains $\gamma_x, \gamma_w \geq 0$ and a rate $\rho \in (0, 1)$ such that 
```{math}
\begin{align*}
\norm{x_k - x^*(x_0)}_2^2 \leq \gamma_x \rho^k \norm{x_0 - x^*(x_0)} + \gamma_w \max_{t \in 0, \ldots, k} \norm{\delta w_t}_2^2, & & \forall k \in \N,  x_0.
\end{align*}
```

Input to state stability is a measure of robustness. It implies linear convergence if $\delta w_k = 0$ for all $k$. 

<!-- The algorithm has a finite  $\ell_2$-gain if there exists a $\gamma_x$ such that 
```{math}
\begin{align*}
 \lim_{T \rightarrow \infty} \frac{\sum_{k=0}^T \norm{ x_k - x^*(x_0)}_2^2}{\sum_{k=0}^T \norm{\delta w_k}_2^2}  \leq \gamma_w   & & \text{if  x_0 = x^*(x_0)}.
\end{align*}
```  -->

More performance criteria are discussed in {doc}`Performance <../usage/problem_formulation/specs>`.

## Well-Posedness

The interconnection $(F, G)$ may be condensed using the set-valued map $H := (F^{-1} - \Dcl)^{-1}$ as
```{math}
\begin{align*}
 (F, G): & & x_{k+1} &\in \Acl x_k + \Bcl H(\Ccl x_k).
\end{align*}
``` 
\
The interconnection is well-posed if  $H$ is globally defined and continuous. Well-posedness implies that the sequence $(x_k, w_k, z_k)_{k \in \N}$ exists and is unique for every  initial state $x_0$. A well-posed algorithm therefore satisfies the explicit relation $x_{k+1} = \Acl x_k + \Bcl H(\Ccl x_k)$ for all $k$. 


{{osyn}} only studies well-posed algorithms. Non-well-posed algorithms such as subgradient descent and Frank-Wolfe are currently out of scope.


## Information Structure and Evaluation

The well-posed algorithmic interconnection $x_{k+1} = \Acl x_k + \Bcl H(\Ccl x_k)$ is a nonlinear iterative procedure. It  may be computationally impractical to compute $x_{k+1}$ from $x_k$ for general $(F, \Dcl)$. The procedure is tractable if $\Dcl$ is block-lower-triangular. The system $(\Acl, \Bcl, \Ccl, \Dcl)$ can be then partitioned as 
```{math}
\begin{align}
\mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2 \\ \vdots \\ z_k^s} &= \mat{c|cccc}{
   \Acl & \Bcl_1 & \Bcl_2  & \hdots & \Bcl_s \hl
   \Ccl_1 & \Dcl_{11} & 0 &   \hdots & 0\\
   \Ccl_2 & \Dcl_{21 } & \Dcl_{22} &   \hdots & 0\\
   \vdots & \vdots & \vdots  & \ddots & \vdots\\
   \Ccl_s & \Dcl_{s1} & \Dcl_{s2} & \hdots & \Dcl_{ss}} \mat{c}{x_{k} \hl w_k^1 \\ w_k^2 \\ \vdots \\ w_k^s}
\end{align}
```


Algorithm evaluation can then causally  proceed for each $k \in \N$ as
```{math}
\begin{align}
w^i_k &= (F_i^{-1} - \Dcl_{ii})^{-1} (\Ccl_i x_k + \textstyle \sum_{j=1}^{i-1} \Dcl_{ij} w^j_k),  & & \forall i \in 1, \ldots, k,\\ 
x_{k+1} &= \Acl x_k + \textstyle \sum_{i=1}^s \Bcl_i w^i_k.
\end{align}
```



The information structure of the algorithm is the block-sparsity pattern of $\Dcl$. If $\Dcl_{ii} = 0$, then $w^i_k$ is explicitly computed from $(x_k, w^1_k, \ldots, w^{i-1}_k)$. If $\Dcl_{ii} \neq 0$, then $w^i_k$ implicitly depends on $(x_k, w^1_k, \ldots, w^{i-1}_k, w^i_k)$. This implicit loop can be broken by proximal/resolvent evaluations for the operator $F_i$, if available. If $\Dcl_{ij} = 0$ with $i > j$, then $w^i_k$ does not use information from the previously computed output $w^j_k$.


Examples of information structures for $s=2$ operators (with $\bullet$ marking  nonzero entries) are 
|    | Sequential    | Parallel       |
| ------: | :-----: | :-------:  |
| Only Implicit | $\mat{cc}{\bullet & 0 \\ \bullet & \bullet}$ | $\mat{cc}{\bullet & 0 \\ 0 & \bullet}$ |
| Mixed | $\mat{cc}{0 & 0 \\ \bullet & \bullet}, \quad  \mat{cc}{\bullet & 0 \\ \bullet & 0}$ | $ \mat{cc}{0& 0 \\ 0 & \bullet}, \quad   \mat{cc}{\bullet & 0 \\ 0 & 0}$ |
| Only Explicit | $\mat{cc}{0 & 0 \\ \bullet & 0}$  | $\mat{cc}{0 & 0 \\ 0 & 0}$ |

The sequential schemes each have a $\bullet$ in the lower-left position: $w^2$ is computed based on information from $w^1$. Parallel schemes can evaluate $w^1$ and $w^2$ separately. 

In analysis, the information structure can be verified by inspection. Synthesis may be constrained to return algorithms with a desired information structure.