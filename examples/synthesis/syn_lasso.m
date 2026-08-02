%% describe the operators
%define the quadratic
rng(430, 'twister');

d = 14;
A_data = randn(2*d, d);
b_data = randn(2*d, 1);

eK = svd(A_data);
m_orig = min(eK)^2;
L_orig = max(eK)^2;

A_norm = A_data/sqrt(m_orig);
b_norm = b_data/sqrt(m_orig);

m = 1;
L = L_orig/m_orig;

op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};



%% form the system
sys = opt_system(ops);

%% solve the problem
config =opt_config();
config.syn.elimination = false;
config.syn.D_mask = [1, 0; 1, 1];
man = opt_synthesis(sys, config);


spec = spec_stability(0.98);
sol = man.solve_single({}, spec);

%% simulate and plot

rng(32, 'twister');

%form the operators
Q = A_norm' * A_norm;
zstar = A_norm \ b_norm;

% Q = rand_quad(d, m, L);
% zstar = 100*(2*rand(d, 1) - 1);
op1_sim = op_sim_quad(Q, zstar);


%define the L1 ball
% tau = 15;
tau = 100;
% op2_sim = op_sim_box(tau);
op2_sim = op_sim_l1_hard(tau);
ops_sim = {op1_sim, op2_sim};

sys_sim = sol.sys.export_sim(ops_sim);
sim = alg_sim(sys_sim, d);
T = 300;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_6f(1);

betastar = sim_out.z(end, :, end);

%% regulator equation tracking


% 
% sim_out_long= sim.sim(3*T); %get a more accurate solution to judge tracking
% 
% betastar = sim_out_long.z(end, :, end);
% wend = sim_out_long.w(:, :, end);
% dstar = [-betastar; wend(1:end-1, :, end)]; %the tracked reference
% 

%find tracking error
% plt = plt.add_opt_sig(sol.regcl, dstar);
% plt.plot_4_err(2);    %plot the tracking error
% plt.plot_4_sq_err(3); %plot the  squared norm of the tracking error