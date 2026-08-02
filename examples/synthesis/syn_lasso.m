%% describe the operators
rng(430, 'twister');

d = 100;
A_data = randn(1.3*d, d);
b_data = 10*rand(1.3*d, 1) + 4*randn(1.3*d, 1);

eK = svd(A_data);
m = min(eK)^2;
L = max(eK)^2;

op1 = op_sml(m, L);
op2 = op_pcc();

ops = {op1, op2};

%% form the system
sys = opt_system(ops);

%% solve the problem
config =opt_config();

%relax stringency of numerical tolerances to encourage a solution

%lowered from default
config.tol.spread = 1e-4;
config.tol.input_diss = 1e-5;
config.tol.M = 1e-9;

%raised from default
config.tol.GX_max = 300;   
config.tol.GY_max = 300;

config.syn.D_mask = [0, 0; 1, 1]; %gradient evaluation of lsq
man = opt_synthesis(sys, config);

sol = man.bisect();

%% simulate and plot

rng(32, 'twister');


%form the operators
op1_sim = op_sim_lsq(A_data, b_data);

%define the L1 ball
tau = 50; %l1 ball constraint
op2_sim = op_sim_l1_hard(tau);

ops_sim = {op1_sim, op2_sim};

sys_sim = sol.sys.export_sim(ops_sim);
sim = alg_sim(sys_sim, d);
sim.sampler.x0 = 7*(2*rand(sys_sim.n, d)-1);
T = 60;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_6f(1);

betastar = sim_out.z(end, :, end);


%% other testing
slift = ss_kron_eye(sol.sys.get_alg(), d);
spartial = lft(A_data'* A_data, slift);
pass_partial = -getPassiveIndex(-spartial, 'input');
szero = lft(zeros(d), spartial);