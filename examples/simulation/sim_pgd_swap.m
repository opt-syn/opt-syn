rng(32, 'twister');

%PGD to minimize quadratic under hard l1 ball constraint
%this uses a time-shifted version of the PGD algorithm


d = 100; %dimension of variable beta

%define the quadratic
m = 1; L = 10;
Q = rand_quad(d, m, L);
zstar = 10*randn(d,1);
op1 = op_sim_quad(Q, zstar);


%define the L1 ball
tau = 100;
op2 = op_sim_l1_hard(tau);
ops = {op2, op1};

%PGD algorithm
gamma = 2/(L + m);
% gamma = 1/L;
K = ss(1, [-gamma, -gamma], [1; 1], [-gamma, 0; -gamma, 0],1);

%add the noise
network = bridge_pass_through(2);
[network, iwp] = network.add_oracle_input(2, []);
[network, izp] = network.perf_output_z(2);


%form the system
sys_clean = opt_system(ops, [], K);
sys_noisy = opt_system(ops, network, K);

%simulate and plot
sim_clean = alg_sim(sys_clean, d);
sim_noisy = alg_sim(sys_noisy, d);
T = 50;
wp_cov = 1;
sim_noisy.sampler.wp = @(k, param) wp_cov*randn(1, d);


sim_clean_out= sim_clean.sim(10*T);
sim_noisy_out = sim_noisy.sim(T);
plt_noisy = alg_plotter(sim_noisy_out);
% plt.plot_4(1);
plt_noisy.plot({'wp', 'w', 'res_w', 'x', 'z', 'f'}, 1);


plt_clean = alg_plotter(sim_clean_out);


%% regulator equation tracking
zend = sim_clean_out.z(:, :, end);

betastar = sim_clean_out.z(end, :, end);
wend = sim_clean_out.w(:, :, end);
dstar = [-betastar; wend(1:end-1, :, end)];
reg = regulator_lti(sys_clean);
regcl = reg.check_regulator();
plt = plt_noisy.add_opt_sig(regcl, dstar);
plt.plot_3_err(2);
plt.plot_3_sq_err(3);