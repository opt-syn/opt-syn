%% describe the operators
rng(430, 'twister');

m = 1;
L = 5;

op1 = op_pcc();
op2 = op_sml(m, L);

ops = {op1, op2};

%% noise corruption in network
network = bridge_pass_through(2);
[network, iwp] = network.add_oracle_input(2, []);
[network, izp] = network.perf_output_z(2);
% [network, izp] = network.perf_output_opt();


%% form the system
sys = opt_system(ops, network);


%% solve the problem
config =opt_config();

%relax stringency of numerical tolerances to encourage a solution

config.syn.reduced_order = true;
config.syn.prox= [1, 0]; %gradient evaluation of lsq
config.syn.elimination = false;
config.recovery.blocks = true;
man = opt_synthesis(sys, config);

%specifications
spec_stab = spec_stability(0.9);
GAIN = 1; %upper bound on gain
spec_stoch = spec_h2(GAIN, 1, 1, 1);
 
spec_stoch.target = true;
% specs = {spec2};
specs = {spec_stab, spec_stoch};
sol = man.solve_single([], specs)
% sol = man.bisect([], specs);

%% simulate and plot
rng(32, 'twister');
d = 40;

%define the L1 ball
tau = 50; %l1 ball constraint
op1_sim = op_sim_l1_hard(tau);

%form the quadratic
Q = rand_quad(d, m, L);
bstar = randi(101, [d, 1]) - 50;
op2_sim = op_sim_quad(Q, bstar);
ops_sim = {op1_sim, op2_sim};

%simulate

sys_sim = sol.sys.export_sim(ops_sim);
sim = alg_sim(sys_sim, d);
sim.sampler.x0 = 7*(2*rand(sys_sim.n, d)-1);
wp_cov = 1;
sim.sampler.wp = @(k, param) wp_cov*randn(1, d);


T = 300;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_6f(1);

sim_clean = sim;
sim_clean.sampler.wp = @(k, param) 0*randn(1, d);
sim_clean_out= sim.sim(T);

% %% other testing
% slift = ss_kron_eye(sol.sys.get_alg(), d);
% spartial = lft(A_data'* A_data, slift);
% pass_partial = -getPassiveIndex(-spartial, 'input');
% szero = lft(zeros(d), spartial);