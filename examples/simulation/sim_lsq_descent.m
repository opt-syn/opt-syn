rng(32, 'twister');

% d = 400; %dimension of variable beta
d = 10;  %number of features
n = 20; %number of constraints

%define the least squares cost
A = randn(n, d);
b = randn(n, 1);
op1 = op_sim_lsq(A, b);

K = (A' * A);
m = min(eig(K));
L = max(eig(K));

zstar = A \ b;
% op1 = op_sim_quad(K, zstar);

ops= {op1};

%PGD algorithm
% gamma = 1/(L + m);
gamma = 2/(L + m);
K = ss(1, [-gamma], [1; ], [0],1);

%form the system
sys = opt_system(ops, [], K);

%simulate and plot
sim = alg_sim(sys, d);
T = 300;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_6f(1);



% %% regulator equation tracking
zend = sim_out.z(:, :, end);
% 
% betastar = sim_out.z(end, :, end);
% wend = sim_out.w(:, :, end);
% dstar = [-betastar; wend(1:end-1, :, end)];
% reg = regulator_lti(sys);
% regcl = reg.check_regulator();
% plt = plt.add_opt_sig(regcl, dstar);
% plt.plot_3_err(2);


