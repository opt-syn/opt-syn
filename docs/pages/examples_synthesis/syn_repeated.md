# Repeated Operators


This example considers composite optimization problems with two functions
\begin{align*}
\beta^* \in \argmin f(\beta) + g(\beta).
\end{align*}

Gradients  $\nabla f$ are 'easy' to evaluate, while subgradients $\partial g$ are computationally expensive and 'hard'. The easy operation $\nabla f$ can be computed $h$ times at every iteration $k$. The hard operation $\partial g$ is computed only once per iteration $k$, and $\partial g$ does not use information from $\nabla f_1$. 

These requirement are formulated using {doc}`information structures <../problem_formulation/config>` and  {doc}`repeated operator evaluation <../usage/problem_formulation/system/bind>`. The associated inclusions and  information structure at $h=3$ is 
\begin{align*}
\mat{c}{w^1_k \\ w^2_k \\ w^3_k \\ w^4_k }  &\in \mat{c}{\nabla f(w^1_k) \\ \nabla f(w^2_k) \\ \nabla f(w^3_k) \\ \partial g(w^4_k)}, & \text{Sparsity}(\Dcl): & \mat{ccc:c}{0 & 0 & 0 & 0 \\ 
\bullet & 0 & 0 & 0 \\
\bullet & \bullet & 0 & 0 \hdl
0 & 0 & 0 & \bullet
}.
\end{align*}

Synthesis is used to find a algorithm that is convergent for all $f \in S_{1, 8}$ and $g \in S_{0,\infty}$, that also satisfies the information structure and repetition requirement. Two rounds of Synthesis/Analysis alternation are performed using `order = {1, 1}`. The Analysis-certified convergence rates for Repeated Evaluations are
:::{list-table}
* - \# $f$ evaluations
  - $\rho$ bound (Round 1)
  - $\rho$ bound (Round 2)  
* - $h=1$
  - 0.7979
  - 0.7982
* - $h=3$
  - 0.9967
  - 0.7050
:::


Figure  [1](#rep) compares convergence behavior of trajectories solving a constrained optimization problem, each starting at $x_0 = 0$.

:::{figure} _static/repeated_compare_dark.png
:align: center
:class: only-dark
:name: rep
*Figure 1:* Convergence conditions for 1-step and 3-step algorithms
:::

:::{figure} _static/repeated_compare_light.png
:align: center
:class: only-light
:name: rep
*Figure 1:* Convergence conditions for 1-step and 3-step algorithms
:::


```{literalinclude} ../../../examples/synthesis/syn_repeated_multigrad.m
:caption: Code for Repeated Synthesis
:language: matlab
:lines:  1-39
```

