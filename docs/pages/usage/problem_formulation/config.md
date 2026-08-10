# Configuration

The configuration options are called by 
``` matlab
config = opt_config();
```

Configuration options include numerical tolerances and  recovery prereferences. These are detailed in the  {doc}`Configuration Documentation <../../documentation/doc_config>`.


## Restricting  Information Structures

The primary user-facing configuration option in algorithm synthesis is specification of the sparsity pattern of the controller matrix $D_K$. 


### Background of Information Structure

The well-posed algorithmic interconnection $x_{k+1} = \Acl x_k + \Bcl H(\Ccl x_k)$ is a nonlinear iterative procedure. Computation of  $x_{k+1}$ from $x_k$ for general $(F, \Dcl)$ requires the solution of a nonlinear fixed-point equation, and may be computationally intractable. If $\Dcl$ is block-lower-triangular, the system $(\Acl, \Bcl, \Ccl, \Dcl)$ can be then partitioned as 
```{math}
\begin{align}
\mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2 \\ \vdots \\ z_k^s} &= \mat{c|cccc}{
   \Acl & \Bcl_1 & \Bcl_2  & \hdots & \Bcl_s \hl
   \Ccl_1 & \Dcl_{11} & 0 &   \hdots & 0\\
   \Ccl_2 & \Dcl_{21 } & \Dcl_{22} &   \hdots & 0\\
   \vdots & \vdots & \vdots  & \ddots & \vdots\\
   \Ccl_s & \Dcl_{s1} & \Dcl_{s2} & \hdots & \Dcl_{ss}} \mat{c}{x_{k} \hl w_k^1 \\ w_k^2 \\ \vdots \\ w_k^s}.
\end{align}
```



The information structure of the algorithm is the block-sparsity pattern of $\Dcl$. 

- If $\Dcl_{ii} = 0$, then $w^i_k$ is explicitly computed from $(x_k, w^1_k, \ldots, w^{i-1}_k)$. 

- If $\Dcl_{ii} \neq 0$, then $w^i_k$ implicitly depends on $(x_k, w^1_k, \ldots, w^{i-1}_k, w^i_k)$.   

- If $\Dcl_{ij} = 0$ with $i > j$, then $w^i_k$ does not use information from the previously computed output $w^j_k$.


Examples of information structures for $s=2$ operators (with $\bullet$ marking  nonzero entries) are 
|    | Sequential    | Parallel       |
| ------: | :-----: | :-------:  |
| Only Implicit | $\mat{cc}{\bullet & 0 \\ \bullet & \bullet}$ | $\mat{cc}{\bullet & 0 \\ 0 & \bullet}$ |
| Mixed | $\mat{cc}{0 & 0 \\ \bullet & \bullet}, \quad  \mat{cc}{\bullet & 0 \\ \bullet & 0}$ | $ \mat{cc}{0& 0 \\ 0 & \bullet}, \quad   \mat{cc}{\bullet & 0 \\ 0 & 0}$ |
| Only Explicit | $\mat{cc}{0 & 0 \\ \bullet & 0}$  | $\mat{cc}{0 & 0 \\ 0 & 0}$ |

The sequential schemes each have a $\bullet$ in the lower-left position: $w^2$ is computed based on information from $w^1$. Parallel schemes can evaluate $w^1$ and $w^2$ separately. 

In Analysis, the information structure can be verified by inspection. Synthesis may be constrained to return algorithms with a desired information structure.

### Control of Information Structures

Coarse-grained control on the sparsity of $K$ is performed by setting the vector `config.syn.prox`, indicating which oracles are permitted proximal evaluation. Fine-grained control is achieved by assigning the matrix `config.syn.D_mask`. 


As an example, a two-operator splitting algorithm with $\Dcl_{11}=0$ (explicit evaluation of $F_1$) and $\Dcl_{22}\neq 0$ (implicit evaluation of $F_2$) can be synthesized using the configuration options
```matlab
%F2 can use information from F1
config.syn.prox = [0, 1];  %method 1 OR
config.syn.D_mask = [0, 0; %method 2
                     1, 1]; 

%F2 cannot use information from F1
config.syn.D_mask = [0, 0; 
                     0, 1]; 
```

Algorithms with both explicit evaluations ($\Dcl_{11}, \Dcl_{22} =  0$) can be synthesized using 
 ```matlab
 %F2 can use information from F1
 config.syn.prox = [0, 0]; %method 1 OR
config.syn.D_mask = [0, 0; %method 2
                      1, 0];

%F2 cannot use information from F1 
config.syn.D_mask = [0, 0; 
                      0, 0]; 
```

By default, `config.syn.D_mask` will be an lower-triangular matrix with all ones, permitting all implicit evaluations. `config.syn.prox` will override `config.syn.D_mask` if both are present.`


An algorithm with nonzero upper-block-triangular  entries of `D_mask` can be Analyzed or Synthesized. However {{osyn}} cannot guarantee that the resulting algorithm will be well-posed, nor will it be able to {doc}`Simulate <../simulation>` trajectories of an algorithm execution.


## Simplified Synthesis Programs

Special structures of the IQC synthesis programs allow for simplification of the LMI programs. 
The supported pairs of simplification methods and dynamical system types are
:::{list-table}
:header-rows: 1
:stub-columns: 1
* -
  - LTI 
  - Periodic-Orbit
  - Periodic
  - Switched
* - Matrix Elimination
  - ✔
  - ✔
  - 
  - 
* - Reduced-Order
  - ✔
  - ✔
  - 
  - 
:::



All simplifications  are enabled by default. They may respectively be disabled by 
```matlab
config.syn.elimination = false;
config.syn.reduced_order = false;
```

### Matrix Elimination

If only one performance requirement is present and the matrix $D_K$ is constrained by `config.syn.D_mask` to have block-lower-triangular sparsity, then a triangular {footcite}`scherer1995complete` Matrix Elimination Lemma {footcite}`gahinet1994linear` can be used to remove the some or all of controller variables from the Synthesis problem. 
A single LMI constraint in Synthesis with the controller variables is replaced by multiple smaller LMI constraints lacking these variables. Both constraint sets have the same feasibility region.

All Synthesis programs the nonlinear  Transformation approach from {footcite}`scherer1997multiobjective` to convexify the search for controllers.

### Reduced-Order Control

Reduced-Order Control uses the internal model structure of the controllers to lower the number of states in the generated algorithm. This reduced-order control is based on the formulation of {footcite}`korouglu2009generalized`. Reduced-order control also allows for the search over solutions of the Regulator Equations.