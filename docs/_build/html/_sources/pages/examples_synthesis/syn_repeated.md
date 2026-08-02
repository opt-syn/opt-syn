# Repeated Operators


This example considers composite optimization problems with two operators
\begin{align*}
\beta^* \in \argmin f(\beta) + g(\beta).
\end{align*}

Gradients  $\nabla f$ are 'easy' to evaluate, while subgradients $\partial g$ are computationally expensive and 'hard'. This requirement is formulated using information structures and  repeated operator evaluation. The easy operation $\nabla f$ can be computed $h$ times at every iteration $k$. The hard operation $\partial g$ is computed only once per iteration $k$, and $\partial g$ does not use information from $\nabla f_1$.


An algorithm allowing for $h=3$ evaluations of $\nabla f_1$ per time index $k$ therefore has the inclusion and information structure
\begin{align*}
\mat{c}{w^1_k \\ w^2_k \\ w^3_k \\ w^4_k }  &\in \mat{c}{\nabla f(w^1_k) \\ \nabla f(w^2_k) \\ \nabla f(w^3_k) \\ \partial g(w^4_k)}, & \text{Sparsity}(\Dcl): & \mat{cccc}{0 & 0 & 0 & 0 \\ 
\bullet & 0 & 0 & 0 \\
\bullet & \bullet & 0 & 0 \\
0 & 0 & 0 & \bullet
}.
\end{align*}

Synthesis is used to find an algorithm satisfying this repetition and information structure requirement. This algorithm must converge for all $f \in S_{1, 8}$ and $g \in S_{0,\infty}$. Two rounds of Synthesis/Analysis alternation are performed using `order = {1, 1}`. Table 

:::{list-table}
:caption: Certified convergence rates by Analysis for Repeated Evaluations
* - \# $f$ evaluations
  - $\rho$ bound (Round 1)
  - $\rho$ bound (Round 2)  
* - $h=1$
  - 0.7979
  - 0.7982
  - $h=3$
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


```{literalinclude} ../../../examples/synthesis/syn_repeated_evaluation.m
:caption: Code for Repeated Synthesis
:language: matlab
:linenos:  1-60
```

