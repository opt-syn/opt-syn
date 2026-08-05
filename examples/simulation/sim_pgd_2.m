rng(32, 'twister');
%https://arxiv.org/pdf/1705.04398
%best rate is 2/(L+m)
%
%but bad transient performance

%PGD to minimize quadratic under hard l1 ball constraint

d = 50; %dimension of variable beta

%define the quadratic
m = 1; L = 100;
Q = rand_quad(d, m, L);
zstar = 100*(2*rand(d, 1) - 1);
op1 = op_sim_quad(Q, zstar);


%define the L1 ball
tau = 200;
op2 = op_sim_l1_hard(tau);
ops = {op1, op2};

%PGD algorithm
gamma = 2/(L + m);
% gamma = 1/L;
K = ss(1, [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);

%form the system
sys = opt_system(ops, [], K);

%simulate and plot
sim = alg_sim(sys, d);
T = 15;
sim.sampler.x0 = 200*randn(1, d);
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
% plt.plot_4(1);
plt.plot({'w', 'res_w', 'z', 'res_z'}, 1);



%% regulator equation tracking
zend = sim_out.z(:, :, end);

betastar = sim_out.z(end, :, end);
wend = sim_out.w(:, :, end);
dstar = [-betastar; wend(1:end-1, :, end)];
reg = regulator_lti(sys);
regcl = reg.check_regulator();
plt = plt.add_opt_sig(regcl, dstar);
plt.plot_3_err(2);