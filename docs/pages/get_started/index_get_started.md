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


A constrained optimization problem minimizing a function $f$ subject to a $L_1$ norm constraint is
```{math}
\begin{align*}
\beta^* \in \text{argmin}_{\norm{\beta}_1 \leq 50} f(\beta).
\end{align*}
```

The function $f$ is known to be real-valued, $m$-strongly convex, and $L$-smooth.

## Analysis 

The Projected Gradient Descent (PGD) algorithm with stepsize $\gamma > 0$ is the iterative procedure

```{math}
\beta_{k+1} = \text{proj}_{\mathcal{Z}}(\beta_k - \gamma \partial f(\beta_k)).
```

PGD achieves linear convergence at rate $\rho \in (0, 1)$ to the optimum $\beta^*$ if there exists a constant $C_0 > 0$ such that
```{math}
\norm{\beta_{k}-\beta^*}_2  \leq  C_0 \rho^{-k} \norm{\beta_0-\beta^*}, \ \forall x_0, k, f.
```

The minimal worst-case convergence rate for PGD is $\frac{L-m}{L+m}$ {footcite}`taylor2018exact`. This rate is achieved with stepsize parameter $\gamma = \frac{2}{L+m}$.


Analysis code to numerically verify this rate for parameters $m = 1, L = 10$ is
```matlab
%describe the operators 
m = 1; L = 10;
rho_theory = (L-m)/(L+m); % 0.8182

op1 = op_sml(m, L); %gradient of f
op2 = op_pcc();     %subdifferential of L1 norm ball indicator function

%State-space representation of PGD
gamma = 2/(L + m);
K = ss(1, [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);

%Interconnect the operators with PGD
sys = opt_system({op1, op2}, [], K);


%run the analysis routine, use bisection to minimize the convergence rate
man = opt_analysis(sys); 
order = {1, 1}; % order of the analysis program
sol = man.bisect(order);
rho = sol.rho % 0.8182, matches PGD theory within 4 digits.
```

A time delay of one time step is introduced before and after evaluation of $\nabla f$. Analysis of PGD with this time delay dynamics with $m=1, L= 10$ is
```matlab
delay = [1, 0];
network = bridge_channel_delay(delay, delay);
sys_delay = opt_system({op1, op2}, network, K);
man_delay = opt_analysis(sys_delay);
sol_delay = man_delay.bisect(order);  % 1.3744
```
Convergence is not guaranteed with time-delays, because $1.3744 > 1$.


## Synthesis

Code to generate an optimization algorithm for $m=1, L=10$ is 
```{literalinclude} ../../examples/getting_started/synthesis_workflow_test.m
:caption: Synthesis without Network Effects
:language: matlab
:lines:  1-13
```

Convergence is confirmed, because the algorithm has a worst-case linear convergence rate of $0.8676 < 1$.

Synthesis is then performed when the oracle $\nabla f$ has a delay of one time step before and after evaluation.
```{literalinclude} ../../../examples/getting_started/synthesis_workflow_test.m
:caption: Synthesis with a 1-step time-delay
:language: matlab 
:lines:  16-21
```

Convergence is again confirmed, since $0.9860 > 1$.

## Simulation

The delay-1 synthesized algorithm is used to solve an $L_1$-norm-constrained quadratic program with $\beta \in \R^{500}$. The code to perform this execution is 
```{literalinclude} ../../examples/getting_started/synthesis_workflow_test.m
:caption: Execution of a 1-step time-delay synthesized algorithm
:language: matlab 
:lines:  25-43
```

Figure [1](#syn-ex) tracing out the execution  is created using the commands
```{literalinclude} ../../../examples/getting_started/synthesis_workflow_test.m
:language: matlab 
:lines:  46-47
```


:::{figure} _static/started_simulation_dark.png
:align: center
:class: only-dark
:name: syn-ex
*Figure 1:* Execution of delay-1 synthesized algorithm
:::

:::{figure} _static/started_simulation_light.png
:align: center
:class: only-light
:name: syn-ex
*Figure 1:*  Execution of delay-1 synthesized algorithm
:::

<!-- The algorithm deployment is convergent to the unique optimum $\beta^* \in \R^{500}$ with 16 nonzero entries. -->

<!-- %0.9860 -->