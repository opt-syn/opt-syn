%% describe the operators
rng(430, 'twister');


m = 1;
L = 10;

op1 = op_sml(m, L);
op2 = op_pcc();

ops = {op1, op2};

%% form the system
sys = opt_system(ops);

%% solve the problem
config =opt_config();

config.syn.D_mask = [0, 0; 1, 1]; %gradient evaluation of lsq
man = opt_synthesis(sys, config);

% sol = man.bisect();
Niter = 4;
order = {[2, 2], [2, 2]};
% order = {[1, 1], [1, 1]};
% order = {3, 3};
% order = {1, 1};
[sol_h, v_h] = man.alternate(Niter, order);

%composite triple-momentum
best_rho = 1 - sqrt(m/L);
[best_rho, v_h(end, end)]

sol = sol_h{1, end}



%% simulate and plot

rng(32, 'twister');


%form the operators
d = 100;
Q = rand_quad(d, m, L);
bstar = 100*randn(d,1);
op1_sim = op_sim_quad(Q, bstar);

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
% slift = ss_kron_eye(sol.sys.get_alg(), d);
% spartial = lft(A_data'* A_data, slift);
% pass_partial = -getPassiveIndex(-spartial, 'input');
% szero = lft(zeros(d), spartial);