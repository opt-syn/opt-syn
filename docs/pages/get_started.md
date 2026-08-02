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
```{literalinclude} ../examples/getting_started/synthesis_workflow_test.m
:caption: Synthesis without Network Effects
:language: matlab
:linenos:  1-13
```

Convergence is confirmed, because the algorithm has a worst-case linear convergence rate of $0.8676 < 1$.



```{literalinclude} ../examples/getting_started/synthesis_workflow_test.m
:caption: Synthesis with a 1-step time-delay
:language: matlab
:linenos:  1-13
```





:::{danger}
This PGD content must be replaced with something else.
:::




## Simulation

We use PGD to solve an optimization problem with $d=400$ dimensions, an $L_1$ norm constraint of $\tau=200$, and a quadratic $f(\beta) = 0.5 (\beta-b_0)'Q (\beta-b_0)$ with $m=1$ and $L = 1000$. 
 The code to execute PGD on this problem for $T=15$ time steps is
```matlab
d = 400; %dimension of variable beta

%define the quadratic (f, randomly generated) 
m = 1; L = 1000;
Q = rand_quad(d, m, L);
b0 = 100*(2*rand(d, 1) - 1);
op1 = op_sim_quad(Q, b0); 


%define the indicator function for the L1 ball (Z)
tau = 200;
op2 = op_sim_l1_hard(tau);

ops = {op1, op2}; %all operators of the problem


%PGD algorithm
gamma = 2/(L + m);
K = ss(1, [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);

%form the system
sys = opt_system(ops, [], K);

%simulate T time steps
sim = alg_sim(sys, d);
T = 15;
sim.sampler.x0 = 200*randn(1, d); %random initial condition
sim_out= sim.sim(T);
```

We then plot the output
```matlab
plt = alg_plotter(sim_out);
plt.plot({'w', 'res_w', 'z', 'res_z'}, 1)
```

:::{figure} _static/sim_pgd_dark.png
:align: center
:class: only-dark
:::

:::{figure} _static/sim_pgd_light.png
:align: center
:class: only-light
:::


## Analysis 


The goal of Analysis is to numerically find the least $\rho$ valid for all $(f, \mathcal{Z})$. We aim to match the worst-case rate $\rho = \frac{L-m}{L+m}$.



The analysis code for parameters $m = 1, L = 10$ is:
``` matlab
%describe the operators 
m = 1; L = 10;
rho_theory = (L-m)/(L+m); % 0.8182

op1 = op_sml(m, L); %\partial f1
op2 = op_pcc();     %\partial f2

%The PGD algorithm using the linear system description
gamma = 2/(L + m);
K = ss(1, [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);

%form the algorithm interconnection
sys = opt_system({op1, op2}, [], K);


%run the analysis routine, use bisection to minimize the convergence rate
man = opt_analysis(sys); 
order = {1, 1}; % order of the analysis program
sol = man.bisect(order);
rho = sol.rho % 0.8182, matches PGD theory within 4 digits.
```

## Synthesis 



The Synthesized optimization algorithm with rate $\rho \leq 0.7209$  is 

```{math}
\begin{align*}
 \mat{c}{x_{k+1}^1 \\ x_{k+1}^2  \hl z_k^1 \\ z_k^2} &= \mat{cc|cc}{I & -0.1326 I & -0.4487 I & 0.3161I  \\
 0 & 0.001327 I &  0.0008995 I & 0.9996 I \hl
  
  I  &  -0.3161 I & -0.3161  I & 0 \\
 I & -0.133 I & -0.449 I & -0.316 I }   \mat{c}{x_{k+1}^1 \\ x_{k+1}^2 \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\partial f(z_k^1) \\ \partial \mathbb{I}_{\mathcal{Z}}(z_k^2)}
\end{align*}.
```


We then perform three rounds of Synthesis-Analysis alternation to acquire an algorithm with certified convergence rate $0.5723 < 0.7209 < 0.8121$.
``` matlab
Niter = 3;
[sol_alt, v_r] = man.alternate(Niter,  {1, 1}); %order {1, 1} Analysis problems
rho_alt = sol_alt{end, end}.rho % 0.5723
```


## Synthesis with Network Dynamics

The  oracle $\partial f$ is no longer directly accessible by the algorithm, but instead takes a nonzero number of steps to interface. 
Algorithm design for a 2-step delay on $\partial f$ is accomplished by executing
``` matlab
%describe the operators 
m = 1; L = 10;
op1 = op_sml(m, L); %\partial f1
op2 = op_pcc();     %\partial f2

%add a round-trip time delay to \partial f1 evaluation
DELAY = [2, 0]; %2-step delay for f1, 0-step (no) delay for f2
network_delay = bridge_channel_delay(DELAY, DELAY); %symmetric delay to and from oracles


%form the algorithm interconnection
sys_delay = opt_system({op1, op2}, network_delay);


%run synthesis, use bisection to minimize the convergence rate
man_delay = opt_synthesis(sys_delay); 
sol = man_delay.bisect();
rho = sol.rho % 0.9553
```