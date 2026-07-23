# Get Started


## Installation

{{osyn}} may be downloaded from [github](https://github.com/Jarmill/opt-syn).

It is tested for MATLAB versions  $\geq$ 2024a.



## Workflow

Analysis and Synthesis follow similar workflows:
1. Define the class of {doc}`operators <usage/problem_formulation/system/operators>' in the optimization/inclusion problem.
2. Specify the algorithm (analysis), or the network interfacing the operators (synthesis)
3. Choose the order of the certification (higher order: better bounds, more expensive)
4. Solve the profiling problem
5. Validate the solution, and plot sample trajectories

## Optimization  Example Setup

We provide an example for analysis and synthesis of a composite optimization algorithm with two operators,
```{math}
\beta^* \in \argmin_{\beta \in \R^c} f_1(\beta) + f_2(\beta),
```

An optimal point $\beta$ must satisfy the necessary inclusion condition $0 \in \partial f_1(\beta) + \partial f_2(\beta)$. 


An instance of this type of composite optimization problem is the [LASSO](https://en.wikipedia.org/wiki/Lasso_(statistics)) task in regression

```{math}
\beta \in \argmin_{||\beta||_1 \leq \tau } f_1(\beta).
```

The class of functions $\F$ we consider in this demonstration are
1. $f_1$ is proper, closed, $m$-strongly convex, and has $L$-Lipschitz gradients with parameters $0 < m < L< \infty$; 
2. $f_2$ is the $0/\infty$ indicator function of a closed, nonempty, convex set $\mathcal{Z}$.

As a result, $\beta^*$ exists and is unique for each pair $(f_1, f_2) \in \F$.

## Analysis 

The Projected Gradient Descent (PGD) algorithm with stepsize $\gamma > 0$ is the iterative procedure

```{math}
\beta_{k+1} = \text{proj}_{\mathcal{Z}}(\beta_k - \gamma \nabla f_1(\beta_k)).
```

PGD achieves linear convergence at rate $\rho \in (0, 1)$ if there exists a constant $\gamma_0$ such that $\norm{\beta_{k}-\beta^*}_2  \leq  \gamma_0 \rho^{-k} \norm{\beta_0-\beta^*}$ for all initial points $\beta_0$, functions $(f_1, f_2) \in \F$, and times $k \in \N$.



The goal of Analysis is to numerically find the least $\rho$ valid for all $f \in \F$. We aim to match the theoretical PGD worst-case linear convergence rate of $\rho = \frac{L-m}{L+m}$, attained by the optimal stepsize $\gamma := \frac{2}{m+L}$. 


The iterative procedure for PGD may be equivalently expressed as an interconnection between a linear dynamical system and the oracles $(f_1, f_2)$ as
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma I & -\gamma I \hl I &0 & 0 \\
 I & -\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\nabla f_1(z_k^1) \\ \partial f_2(z_k^2)}
\end{align*},
```

The algorithm is convergent if  $\lim_{k \rightarrow \infty} z_1 = \lim_{k \rightarrow \infty} z_2 = \beta^*$ for all  $(f_1, f_2) \in \F$. 

The analysis code at $m = 1, L = 10$ is:
``` matlab
%describe the operators 
m = 1; L = 10;
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
rho = sol.rho % 0.8182, matches theory within 4 digits.
```

## Synthesis 

In synthesis, we aim to design an algorithm with minimal convergence rate. The synthesis code for the case of $m=1, L=10$ is 

``` matlab
%describe the operators 
m = 1; L = 10;
op1 = op_sml(m, L); %\partial f1
op2 = op_pcc();     %\partial f2

%form the algorithm interconnection
sys = opt_system({op1, op2});


%run the synthesis routine, use bisection to minimize the convergence rate
man = opt_synthesis(sys); 
sol = man.bisect();
rho = sol.rho % 0.7209
```

The generated optimization $\rho \leq 0.7209$  algorithm from synthesis is 

```{math}
\begin{align*}
 \mat{c}{x_{k+1}^1 \\ x_{k+1}^2  \hl z_k^1 \\ z_k^2} &= \mat{cc|cc}{I & -0.1326 I & -0.4487 I & 0.3161I  \\
 0 & 0.001327 I &  0.0008995 I & 0.9996 I \hl
  
  I  &  -0.3161 I & -0.3161  I & 0 \\
 I & -0.133 I & -0.449 I & -0.316 I }   \mat{c}{x_{k+1}^1 \\ x_{k+1}^2 \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\nabla f_1(z_k^1) \\ \partial f_2(z_k^2)}
\end{align*}.
```


We then perform three rounds of Synthesis-Analysis alternation to acquire an algorithm with certified convergence rate $0.5723 < 0.8121$.
``` matlab
b_opts = bisect_opts;
b_opts.Niter = 3;
sol_alt = man_grad.alternate([], {1, 1}, [], b_opts);
rho_alt = sol_alt{end, end}.rho % 0.5723
```

<!-- 
These algorithm require proximal evaluation of the function $f_1$,  due to the nonzero quantity $-0.3161$ between $w_1$ and $z_1$. If only gradient evaluations of $\nabla f_1$ are permitted, then we use the synthesis code

``` matlab
config_grad = opt_config();
config_grad.syn.D_mask = [0, 0; 1, 1]; %enforce only gradient evaluation of f1
man_grad = opt_synthesis(sys, config_grad);
sol_grad = man_grad.bisect();
rho_grad = sol_grad.rho % 0.8322
```

to obtain the $\rho \leq 0.8322$ algorithm 
```{math}
\begin{align*}
 \mat{c}{x_{k+1}^1 \\ x_{k+1}^2  \hl z_k^1 \\ z_k^2} &= \mat{cc|cc}{I & -0.05622 I & -0.1715 I & -0.2277 I \\
 0 &  - 0.0008595 I & -0.00156 I & 1.007 I \hl
  
  I  & 0 &  0 & 0 \\
 I & 0.3658 I & -0.1153 I & -0.4811 I }   \mat{c}{x_{k+1}^1 \\ x_{k+1}^2 \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{\nabla f_1(z_k^1) \\ \partial f_2(z_k^2)}
\end{align*}.
``` -->



## Synthesis with Network Dynamics

The gradient oracle $\nabla f_1$ is no longer directly accessible by the algorithm, but instead takes a nonzero number of steps to interface. 
Algorithm design for a 2-step delay on $\nabla f_1$ is accomplished by executing
``` matlab
%describe the operators 
m = 1; L = 10;
op1 = op_sml(m, L); %\partial f1
op2 = op_pcc();     %\partial f2

%add a round-trip time delay to f1 evaluation
DELAY = [2, 0]; $2-step delay for f1, 0-step (no) delay for f2
network_delay = bridge_channel_delay(DELAY, DELAY); %symmetric delay to and from oracles


%form the algorithm interconnection
sys_delay = opt_system({op1, op2}, network_delay);


%run the synthesis routine, use bisection to minimize the convergence rate
man_delay = opt_synthesis(sys_delay); 
sol = man_delay.bisect();
rho = sol.rho % 0.9553
```

## Validation

The `sol` structure contains information about the solution of analysis/synthesis. The solution is feasible if the following conditions are met

| Name   |  Description  | Valid Condition |
|----| ---- | ----- | 
| `STATUS` | Feasibility of problem | 0 if feasible, nonzero if infeasible |
| `dia` | Constraint violation | `dia`<0 if strictly feasible, `dia`=0 if marginally feasible, `dia` > 0 if infeasible |
| `gain` | Input passivity index and $H_\infty$ gain | Feasible if `gain(1)` < 0 and `gain(2)` < 1 |

If all of the above conditions are met, then linear convergence is established if and only if `sol.rho` < 1. A finite `sol.rho` > 1 establishes a bounded rate of divergence. No conclusions can be drawn about linear convergence if `sol.rho` = 1. 

<!-- ## Plotting

After these numerical checks, the  -->