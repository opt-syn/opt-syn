#  Channel Memory

This example involves a two-operator inclusion problem. Memory effects are present in the communication link to and from evaluation of $F_1$. 
The intensity of the memory effects are represented by a scalar $\alpha > 0$. 
The network effects are described by a transfer function matrix
```{mathbf}
 \mat{c}{z \hdl y} =  \mat{cc:cc}{0 & 0 & \frac{\mathbf{z}}{\mathbf{z} +\alpha} I& 0 \\
0 & 0 & 0 & I \hdl
\frac{\mathbf{z}}{\mathbf{z} +\alpha} I& 0 & 0 & 0 \\
0 & I & 0 & 0}  \mat{c}{w \hdl u},
```
where $\mathbf{z}$ is the unit shift operator $((\mathbf{z} x)_k :=  x_{k+1}$ for all $k \in \N$).  

A state-space realization of this transfer function matrix is 
```{math}
 \text{Network}: \qquad  \mat{c}{x_{k+1}^N \hl z_k \hdl y_k} =  \mat{cc|cc:cc}{-\alpha I & 0 & \frac{1}{2} I & 0 & 0 & 0 \\
 0 & -\alpha I & 0 & 0 & \frac{1}{2}I & 0  \hl
 0 & -2\alpha I & 0 & 0 & I & 0  \\
 0 & 0 & 0 & 0 & 0 & I  \hdl
 -2 \alpha I & 0 & I & 0 & 0 & 0 \\
 0 & 0 & 0 & I & 0 & 0}  \mat{c}{x_{k}^N \hl w_k \hdl u_k}.
```

The degenerate case of $\alpha=0$ is no network dynamics.

The controller structure with parameters $(\gamma, \lambda) \geq 0$ used to solve the two-operator inclusion problem is
```{math}
 \text{Network}: \qquad  \mat{c}{x_{k+1}^c \hl u_k^1 \\ u_k^2} =  \mat{c|cc}{I & -\gamma \lambda & -\gamma \lambda \frac{1}{\alpha+1} \\
 (1+\alpha) I & 0 & 0 \\
 I & -\gamma I & -\gamma\frac{1}{\alpha+1} I }  \mat{c}{x_{k}^c \hl y_k^1 \\ y_k^2}.
```

This controller structure is parameterized by $\alpha$. If $\alpha = 0$ and $\lambda = 1$, then this controller is the same as Projected Gradient Descent.
The controller structure is chosen to ensure that the {doc}`Regulator Equation <../how_it_works/network_synthesis>` condition for algorithm convergence is satisfied for all values $(\gamma, \lambda). Projected Gradient Descent fails the regulator equation requirement of convergence when $\alpha > 0$. 

We use this algorithm to solve a composite optimization problem 
```{math}
\beta^* \in \argmin_{\norm{\beta}_1 \leq 100 } f(\beta)
```
where $f$ is a convex quadratic with $m=1$, $L=5$. 

Figure [1](#sym-trace) plots a trace of algorithm execution starting from $x_0 = 0$, highlighting the states of the network and controller.

:::{figure} _static/sim_channel_sym_dark.png
:align: center
:class: only-dark
:name: sym-trace
*Figure 1:* trace of execution and convergence
:::

:::{figure} _static/sim_channel_sym_light.png
:align: center
:class: only-light
:name: sym-trace
*Figure 1:* trace of execution and convergence
:::

We now establish tracking properties using the Regulator Equations. The problem has a unique optimal solution $\beta^*$ since $f$ is strongly convex, and the constraint set is nonempty. Because $f$ is smooth, there is a unique vector $w^{*1} = \nabla f(\beta^*)$, and therefore the subdifferential at optimality $w^* = (\nabla f(\beta^*),  -\nabla f(\beta^*))$ is likewise unique. 

The solution to the Regulator Equations for this network and controller are
```{math}
\Pi &= \mat{cc}{0 & \frac{1}{2(\alpha+1)}I \\ -\frac{1}{2} I & 0 }, \Gamma &= \mat{cc}{-(\alpha+1) I & 0 \\
-I & 0 }, & \Phi &= \mat{cc}{0 & \frac{1}{2(\alpha+1)}I \\ 0  & -I }, \Theta &= \mat{c}{-I & 0}.
```

Algorithm convergence implies tracking of the signals $(x^N, x^c, y, u)$ with
```{math}
\begin{align}
    \lim_{k \rightarrow \infty}
    \mat{c}{x_k^N \\ x^c_k \hl  y_k \\u_k} & = \mat{c}{\Pi \\ \Theta \hl \Phi \\ \Gamma} \mat{c}{-\beta^* \\ \nabla f(\beta^*)}.
\end{align}
```


Figure [2](#sym-track) plots the tracking of these signals over time

:::{figure} _static/sim_channel_sym_track_dark.png
:align: center
:class: only-dark
:name: sym-track
*Figure 2:* Tracking error
:::

:::{figure} _static/sim_channel_sym_track_light.png
:align: center
:class: only-light
:name: sym-track
*Figure 2:* Tracking error
:::

Figure [3](#sym-trace-sq) plots the squared norm of the tracking error


:::{figure} _static/sim_channel_sym_track_sq_dark.png
:align: center
:class: only-dark
:name: sym-track-sq
*Figure 3:* Tracking residuals
:::

:::{figure} _static/sim_channel_sym_track_sq_light.png
:align: center
:class: only-light
:name: sym-track-sq
*Figure 3:* Tracking residuals
:::

Code to generate this simulation is
```matlab
rng(32, 'twister');

d = 200; %dimension of variable beta

%% describe the operators
%define the quadratic
m = 1; L = 5;
Q = rand_quad(d, m, L);
zstar = 100*(2*rand(d, 1) - 1);
op1 = op_sim_quad(Q, zstar);


%define the L1 ball
tau = 100;
op2 = op_sim_l1_hard(tau);
ops = {op1, op2};

gamma = 0.4;
lambda = 0.2;

%% describe the network
%use a transfer function representation 
%to model the channel memory
z = tf('z', 1);
alpha = 0.4; %channel memory effect
ascale = (2*(alpha+1));

P = [0, 0, z/(z+alpha), 0;
    0 , 0, 0, 1;
    z/(z+alpha), 0, 0, 0;
    0, 1, 0, 0];

%partition the input and output channels
network = genplant(P);
network.nw = 2; network.nu = 2;
network.nz = 2; network.ny = 2;

%% describe the controller
K = ss([1], [-gamma*lambda, -gamma*lambda/(alpha+1)], ...
    [1+alpha; 1], [0, 0; -gamma, -gamma/(alpha+1)],1);


%% form the system
sys = opt_system(ops, network, K);


%% simulate and plot
sim = alg_sim(sys, d);
T = 100;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_6(1);




%% plot tracking error

%get a more accurate solution to judge tracking error
sim_out_long= sim.sim(3*T); 

%compute the tracked reference 
betastar = sim_out_long.z(end, :, end);
wend = sim_out_long.w(:, :, end);
dstar = [-betastar; wend(1:end-1, :, end)]; 


%compute the regulator equation solution
reg = regulator_lti(sys);      %solution to open-loop regulator equations
regcl = reg.check_regulator(); %solution to closed-loop regulator equations

plt = plt.add_opt_sig(regcl, dstar);
plt.plot_4_err(2);    %plot the tracking error
plt.plot_4_sq_err(3); %plot the  squared norm of the tracking error
```