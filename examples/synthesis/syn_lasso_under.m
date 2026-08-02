% lasso in the underparameterized setting

%% describe the operators
%define the quadratic
rng(430, 'twister');

% d = 10; %dimensionality
d = 300;
tau = 50; %l1 ball constraint
A_data = rand(d/2, d);
b_data = rand(d/2, 1);

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
% config.syn.D_mask = [0, 0; 1, 1]; %gradient evaluation of lsq
config.syn.D_mask = [1, 0; 1, 1];   %prox evaluation of lsq


%for sublinear convergence, set lower bounds to the dissipation terms to
%zero (nonstrict dissipation)
config.tol.input_diss = 0;
config.tol.M = 0;
% config.tol.spread = 1e-3;
man = opt_synthesis(sys, config);
% man = opt_synthesis(sys, config);


sol = man.solve_single({}, spec_erg);

%% simulate and plot


%form the operators
op1_sim = op_sim_lsq(A_data, b_data);

%define the L1 ball

op2_sim = op_sim_l1_hard(tau);

ops_sim = {op1_sim, op2_sim};

sys_sim = sol.sys.export_sim(ops_sim);
sim = alg_sim(sys_sim, d);
T = 1000;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_6f(1);

betastar = sim_out.z(end, :, end);



