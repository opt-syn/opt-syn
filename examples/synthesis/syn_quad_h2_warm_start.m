%% describe the operators
rng(430, 'twister');

m = 1; L = 10;

op1 = op_pcc();
op2 = op_sml(m, L);

ops = {op1, op2};

%% noise corruption in network
network = bridge_pass_through(2);
[network, iwp] = network.add_oracle_input(2, []);
[network, izp] = network.perf_output_z(2);


%% form the system
sys = opt_system(ops, network);


%% solve the problem
config =opt_config();
config.syn.prox= [1, 0]; %prox evaluation of g, gradient evaluation of f
config.syn.elimination = true;
man = opt_synthesis(sys, config);

%specifications
spec_stab = spec_stability(0.9);
GAIN = 10; %upper bound on gain
spec_stoch = spec_h2(GAIN, 1, 1, 1);
 
spec_stoch.target = true;
specs = {spec_stab, spec_stoch};

%create warm start
gamma = 2/(L + m);
K = ss(1, [-gamma, -gamma], [1; 1], [-gamma, 0; -gamma, 0],1);
order = 1;
sys_ana = sys;
sys_ana.K = K;
man_ana = opt_analysis(sys_ana, config);
sol_ana = man_ana.bisect(order); %rho = 0.8182
warm_start = sol_ana.cert.iqc_op;

%only verify linear convergence
sol_exp = man.bisect(); %rho = 0.8321
sol_exp_warm = man.bisect(warm_start); %rho = 0.7476


%optimize over convergence rate (outer) and stochastic sensitivity (inner)
sol = man.bisect([], specs); %rho = 0.8561, h2 = 0.3792
sol_warm = man.bisect(warm_start, specs); %rho = 0.8395, h2 = 0.9737



%% attempt alternation
iqc_warm =  sol_joint.cert.iqc_op; %from ana_pgd_h2;
sol_warm = man.bisect(iqc_warm, specs);
sol_single = man.solve_single(iqc_warm, specs);
% order = {2, 2};
% Niter = 3;
% [sol_alt, v_h]= man.alternate(Niter, order, [], specs);

%% simulate and plot
rng(32, 'twister');
d = 40;

%define the L1 ball
tau = 50; %l1 ball constraint
op1_sim = op_sim_l1_hard(tau);

%form the quadratic
Q = rand_quad(d, m, L);
bstar = randi(101, [d, 1]) - 50;
% bstar = zeros(d, 1);
op2_sim = op_sim_quad(Q, bstar);
ops_sim = {op1_sim, op2_sim};

%simulate
sys_sim = sol.sys.export_sim(ops_sim);
sim = alg_sim(sys_sim, d);
sim.sampler.x0 = 7*(2*rand(sys_sim.n, d)-1);
wp_std = 10;
sim.sampler.wp = @(k, param) wp_std*randn(1, d);


T = 101;
% T = 2000;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_6f(1);
% plt.plot({'wp', 'res_w', 'z', 'f'}, 1);



%% regulator

sys_sim = sol.sys.export_sim(ops_sim);
sim2 = sim;
sim2.sampler.wp = @(k, param) 0*randn(1, d);

sim_clean = sim;
sim_clean.sampler.wp = @(k, param) 0*randn(1, d);
sim_clean_out= sim.sim(T);
sim2_out= sim.sim(T);
betastar = sim2_out.z(end, :, end);
wstar = sim2_out.w(end, :, end);
dstar = [-betastar; wstar];
reg = regulator_lti(sys_sim);
regcl = reg.check_regulator();
plt = plt.add_opt_sig(regcl, dstar);


%% plot the error

zerr2 = sim_out.z(end, :, :) - betastar;
sq_zerr2 = squeeze(sum(zerr2.^2, 2));
figure(5)
clf
hold on
cc = linspecer(2);
plot(sim_out.k, sq_zerr2, 'linewidth', 2, 'color', cc(1, :))
mean_bound = d*wp_std^2 * (sol.objective);
mean_emp = mean(sq_zerr2(200:end));
xl = xlim;

plot(xl, mean_bound*[1, 1], ':', 'color', 0.5*[1,1,1], 'linewidth', 2);
plot(xl, mean_emp*[1, 1],  '--', 'linewidth', 2, 'color', cc(2, :));
set(gca, 'yscale', 'log')
ylim([0.5*min(sq_zerr2), 1.5*mean_bound])


xlabel('$k$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$||z^2 - \beta^*||^2_2$', 'interpreter', 'latex', 'fontsize', 16)
rho_best = (L-m)/(L+m);
xl = xlim;
plot(xl, rho_best * [1, 1], ':', 'color', 0.5 * [1, 1, 1], 'LineWidth', 2);
lname = {'Performance Output', 'Mean (bound)', 'Mean (empirical)'};
legend(lname, 'location','northeast', 'interpreter', ...
    'latex', 'fontsize', 16)


% %% other testing
% slift = ss_kron_eye(sol.sys.get_alg(), d);
% spartial = lft(A_data'* A_data, slift);
% pass_partial = -getPassiveIndex(-spartial, 'input');
% szero = lft(zeros(d), spartial);