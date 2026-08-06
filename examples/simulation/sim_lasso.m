rng(32, 'twister');

d = 400; %dimension of variable beta
% d = 10;  %number of features
n = 15; %number of constraints

%define the least squares cost
A = randn(n, d);
b = randn(n, 1);
op1 = op_sim_lsq(A, b);


Q = (A' * A);
m = min(eig(Q));
L = max(eig(Q));

% op1 = op_sim_quad(Q, A \ b);

%define the L1 ball
% op2 = op_sim_box(10);
tau = 50;
op2 = op_sim_l1_hard(tau);
% op2 = op_sim_quad(Q, ones(size(A \ b)));
ops = {op1, op2};

%PGD algorithm
% gamma = 2/(L + m) + 1e-4;
% a = 0.03;
% a = 1e-4;
% a = 0;
% a = 0.5;
a = 1;
gamma = 1/L * a  + (1-a) * 2/(L+m);
% gamma = 1/L;
K = ss(1, [-gamma, -gamma], [1; 1], ...
    [0, 0; -gamma, -gamma],1);
% K = ss(1, [-gamma, -gamma], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);

%form the system
sys = opt_system(ops, [], K);

%simulate and plot
sim = alg_sim(sys, d);
% T = 60;
% T = 400;
T = 500;
% T = 1000;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_4(1);

% 
% %% take a different path
% 
% gamma = 1/L;
% % gamma = 1/L;
% K2 = ss(1, [-gamma, -gamma], [1; 1], ...
%     [0, 0; -gamma, -gamma],1);
% sim2 = sim
% sim.sys.K = K2;
% 
% sim_out2= sim.sim(500);
% plt2 = alg_plotter(sim_out2);
% plt2.plot_4(2);
% 
% %% regulator equation tracking
% betastar = sim_out2.z(end, :, end);
% wend = sim_out.w(:, :, end);
% dstar = [-betastar; wend(1:end-1, :, end)];
% reg = regulator_lti(sys);
% regcl = reg.check_regulator();
% plt2 = plt2.add_opt_sig(regcl, dstar);
% plt2.plot_3_err(3);
% 
% %% bind to the first one
% 
% plt = plt.add_opt_sig(regcl, dstar);
% plt.plot_3_err(3);