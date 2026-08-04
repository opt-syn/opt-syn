% lasso in the underparameterized setting

%% describe the operators
%define the data matrix
rng(430, 'twister');
d = 200;
A_data = 0.5*randn(d/2, d);
b_data = 10*randn(d/2, 1) + 2*rand(d/2, 1);

eK = svd(A_data);
m = 0;
L = max(eK)^2;

op1 = op_sml(m, L);
op2 = op_sml(m, L);
ops = {op1, op2};



%% form the system
sys = opt_system(ops);
spec_erg = spec_ergodic();

%% solve the problem
config =opt_config();
config.syn.D_mask = [0, 0; 1, 1]; %gradient evaluation of lsq

%for sublinear convergence, set lower bounds to the dissipation terms to
%zero (nonstrict dissipation)
config.tol.input_diss = 0;
config.tol.M = 0;
man = opt_synthesis(sys, config);

sol = man.solve_single({}, spec_erg);

%% simulate and plot


%form the operators
op1_sim = op_sim_lsq(A_data, b_data);

%define the L1 ball
tau = 50; %l1 ball constraint
op2_sim = op_sim_l1_hard(tau);

ops_sim = {op1_sim, op2_sim};

sys_sim = sol.sys.export_sim(ops_sim);
sim = alg_sim(sys_sim, d);
sim.sampler.x0 = 7*(2*rand(sys_sim.n, d)-1);
T = 10000;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_6f(1);

betastar = sim_out.z(end, :, end);

betaols = A_data \ b_data;