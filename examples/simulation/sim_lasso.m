rng(32, 'twister');

% d = 400; %dimension of variable beta
d = 10;  %number of features
n = 15; %number of constraints

%define the least squares cost
A = randn(n, d);
b = randn(n, 1);
op1 = op_sim_lsq(A, b);

K = (A' * A);
m = min(eig(K));
L = max(eig(K));


%define the L1 ball
% op2 = op_sim_box(40);
tau = 100;
op2 = op_sim_l1_hard(tau);
ops = {op1, op2};

%PGD algorithm
% gamma = 2/(L + m);
gamma = 1/L;
K = ss(1, [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);
% K = ss(1, [-gamma, -gamma], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);

%form the system
sys = opt_system(ops, [], K);

%simulate and plot
sim = alg_sim(sys, d);
T = 100;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_4(1);



% %% regulator equation tracking
% zend = sim_out.z(:, :, end);
% 
% betastar = sim_out.z(end, :, end);
% wend = sim_out.w(:, :, end);
% dstar = [-betastar; wend(1:end-1, :, end)];
% reg = regulator_lti(sys);
% regcl = reg.check_regulator();
% plt = plt.add_opt_sig(regcl, dstar);
% plt.plot_3_err(2);