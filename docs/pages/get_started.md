# Get Started


## Installation

{{osyn}} may be downloaded from [github](https://github.com/Jarmill/opt-syn).

It is tested for MATLAB versions  $\geq$ 2024a.



## Workflow

Analysis and Synthesis follow similar workflows:
1. Define the class of functions/operators in the optimization/inclusion problem.
2. Specify the algorithm (analysis), or the network interfacing the operators (synthesis)
3. Choose the order of the certification (higher order: better bounds, more expensive)
4. Solve the profiling problem
5. Validate the solution, and plot sample trajectories

(#optimization-example-setup)=
## Optimization  Example Setup


A constrained optimization problem of minimizing a function $f$ subject to a sparse $L_1$ norm constraint is
```{math}
\begin{align}
\beta^* \in \text{argmin}_{\norm{\beta}_1 \leq 100} f(\beta), & & & \beta^* \in \text{argmin}_{\beta \in \R^d} f(\beta) + \mathbf{I}_{\norm{\cdot}_1 \leq 100}(\beta),
\end{align}
```

The function $f$ is known to be real-valued, $1$-strongly convex, and $50$-smooth (Lipschitz gradients).


$f$ only accessible by the optimizer using a network, which may possess time delays or other dynamics.  The goal is to Synthesize  certifiably convergent algorithms that will solve the optimization problem in this remote environment. 




## Synthesis

When the function $f$ is directly connected to the optimizer (no network dynamics), the code to perform synthesis is
```{literalinclude} ../../examples/getting_started/synthesis_workflow_test.m
:caption: Synthesis without Network Effects
:language: matlab
:lines:  1-13
```

Convergence is confirmed, because the algorithm has a worst-case linear convergence rate of $0.8676 < 1$.

Synthesis is then performed when the oracle $\nabla f$ has a time delay of one step in each direction.

```{literalinclude} ../../examples/getting_started/synthesis_workflow_test.m
:caption: Synthesis with a 1-step time-delay
:language: matlab %0.9860
:lines:  16-21
```

:::{caution} 
Under Construction
:::