# Noisy Douglas Rachford


The Douglas-Rachford algorithm is a procedure for solving a two-operator inclusion problem {footcite}`douglas1956numerical`.
It is characterized by  parameters $\gamma, \lambda > 0$, and can be described as the interconnection
```{math}
\begin{align*}
 \mat{c}{x_{k+1} \hl z_k^1 \\ z_k^2} &= \mat{c|cc}{I & -\gamma \lambda I & -\gamma \lambda I \hl I &-\gamma I & 0 \\
 I & -2\gamma I & -\gamma I  }   \mat{c}{x_{k} \hl w_k^1 \\ w_k^2}, & \mat{c}{w_k^1 \\ w_k^2} \in  \mat{c}{F_1(z_k^1) \\  F_2(z_k^2)}.
\end{align*}
```



We  use Douglas-Rachford to solve a composite optimization problem
```{math}
\beta^* \in \argmin_{\norm{\beta}_\infty \leq 10} f(\beta)
```

We describe the necessary optimality condition for this problem using the operators $F_1 = \partial f$, and $F_2 = \partial \mathbb{I}_{\norm{\cdot}_\infty \leq 10}.$ 


## No Noise
We execute the Douglas-Rachford scheme with parameters $\gamma = 0.4, \lambda = 1$ for a problem where $f_1$ is a convex quadratic (eigenvalue bounds $m = 1, L = 10$). Figure [1](#dr-clean) plots a trajectory of Douglas-Rachford starting from $x_0 = 0$ is

:::{figure} _static/dr_clean_dark.png
:align: center
:class: only-dark
:name: dr-clean
:::

:::{figure} _static/dr_clean_light.png
:align: center
:class: only-light
:name: dr-clean
:::


## With Noise
We then add noise to the Douglas-Rachford scheme. The performance input $w_p$ is additive noise
at the output of the subgradient evaluations. The performance output $z_p$ is the consensus error $\frac{1}{2}(z_1 - z_2)$ at each channel. 
The {doc}`System <../usage/problem_formulation/system/index_system>` representing Douglas-Rachford with this noise structure is
```{math}
\begin{align*}
 \text{Operator} & & \mat{c}{w_k^1 \\ w_k^2} &\in  \mat{c}{F_1(z_k^1) \\  F_2(z_k^2)}, \\
\text{Network} & & \mat{c}{z_k^1 \\ z_k^2 \hdl z_{p, k}^1 \\ z_{p, k}^2 \hdl y_{k}^1\\y_k^2} &= \mat{cc:cc:cc}{0 & 0 & 0 & 0 & I & 0 \\
0 & 0 & 0 & 0 & 0 & I  \hdl
0 & 0 & 0 & 0 & \frac{1}{2} I & -\frac{1}{2} I \\
0 & 0 & 0 & 0 & -\frac{1}{2} I & \frac{1}{2} I \hdl
0 & 0 & I & 0 & I & 0 \\
0 & 0 & 0 & I & 0 & I} \mat{c}{w_k^1 \\ w_k^2 \hdl w_{p, k}^1 \\ w_{p, k}^2 \hdl u_{k}^1 \\ u_k^2}, \\
 \text{Douglas-Rachford} & & \mat{c}{x^c_{k+1} \hl u_k^1 \\ u_k^2} &= \mat{c|cc}{I & -\gamma \lambda I & -\gamma \lambda I \hl I &-\gamma I & 0 \\
 I & -2\gamma I & -\gamma I  }   \mat{c}{x^c_{k} \hl y_k^1 \\ y_k^2}.
\end{align*}
```

Figure [2](#dr-noisy) plots a trace of a trajectory starting at $x_0$, in which the performance input $w_p$ is randomly sampled subject to the bound $\norm{w_{p,k}}_2 \leq 10, \forall k \in \N$.
:::{figure} _static/dr_noisy_dark.png
:align: center
:class: only-dark
:name: dr-noisy
Response under bounded noise
:::

:::{figure} _static/dr_noisy_light.png
:align: center
:class: only-light
:name: dr-noisy
Response under bounded noise
:::


## Code

The code to generate this demonstration is 
```matlab
%Douglas Rachford Algorithm with noise
rng(32, 'twister');

d = 100; %dimension of variable beta
s = 2; %number of operators

%% create the system
%define the quadratic
m = 1; L = 10;
Q = rand_quad(d, m, L);
zstar = 100*(2*rand(d, 1) - 1);
op1 = op_sim_quad(Q, zstar);

%define the L infinity ball
BOX = 10;
op2 = op_sim_box(BOX);
ops = {op1, op2};

%douglas-rachford
gamma = 0.4;  lambda = 1;   %stepsizes
K = ss(1, [-gamma*lambda, -gamma*lambda], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);



%system with no noise
sys_clean = opt_system(ops, [], K);




%now create the network for the system with noise
network = bridge_pass_through(s);
%add noise to the subgradients (outputs of oracles)
iwp = [1,2]; izp = [];
network = network.add_oracle_input(iwp, izp);

%consensus error as a performance condition
network  = network.perf_output_con();

%form the system with noise
sys = opt_system(ops, network, K);

%% simulate and plot
T = 100; %time horizon

%no noise
sim_clean = alg_sim(sys_clean, d);
sim_out_clean = sim_clean.sim(T);
plt_clean = alg_plotter(sim_out_clean);
fig_clean = plt_clean.plot({'x', 'w', 'res_w', 'f', 'z', 'res_z'},  2);

%with noise
sim = alg_sim(sys, d);

%perform sampling
eps_w = 10; %norm(w_k, 2) <= eps_w at each k
sim.sampler.wp = @(param)  eps_w * (ball_sample(length(iwp), d));

sim_out = sim.sim(T);
plt_noisy = alg_plotter(sim_out);
fig_noisy = plt_noisy.plot({'wp','w', 'res_w', 'zp','z',  'res_z'},  3);
```


