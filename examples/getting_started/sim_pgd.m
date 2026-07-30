rng(32, 'twister');

%PGD to minimize quadratic under hard l1 ball constraint

d = 400; %dimension of variable beta

%define the quadratic
m = 1; L = 1000;
Q = rand_quad(d, m, L);
zstar = 100*(2*rand(d, 1) - 1);
op1 = op_sim_quad(Q, zstar);


%define the L1 ball
tau = 200;
op2 = op_sim_l1_hard(tau);
ops = {op1, op2};

%PGD algorithm
gamma = 2/(L + m);
K = ss(1, [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);

%form the system
sys = opt_system(ops, [], K);

%simulate and plot
sim = alg_sim(sys, d);
T = 15;
sim.sampler.x0 = 200*randn(1, d);
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
% plt.plot_4(1);
plt.plot({'w', 'res_w', 'z', 'res_z'}, 1);