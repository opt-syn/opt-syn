rng(32, 'twister');

%Douglas Rachford Algorithm with noise

d = 100; %dimension of variable beta
s = 2; %number of operators

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

%create the network
network = bridge_pass_through(s);
%add noise to the subgradients (outputs of oracles)
iwp = [1,2];
izp = [];
network = network.add_oracle_input(iwp, izp);

%consensus error as a performance condition
network  = network.perf_output_con();

%form the system
sys = opt_system(ops, network, K);
sys_clean = opt_system(ops, [], K);

%% simulate and plot
T = 100;

%no noise
sim_clean = alg_sim(sys_clean, d);
sim_out_clean = sim_clean.sim(T);
plt_clean = alg_plotter(sim_out_clean);
fig_clean = plt_clean.plot({'x', 'w', 'res_w', 'f', 'z', 'res_z'},  2);

%with noise
sim = alg_sim(sys, d);

eps_w = 10; %norm(w_k, 2) <= eps_w at each k
sim.sampler.wp = @(param)  eps_w * (ball_sample(length(iwp), d));

sim_out = sim.sim(T);
plt_noisy = alg_plotter(sim_out);
fig_noisy = plt_noisy.plot({'wp','w', 'res_w', 'zp','z',  'res_z'},  3);