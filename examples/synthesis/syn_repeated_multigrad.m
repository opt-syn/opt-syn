rng(33, 'twister');

m = 1;
L = 8;
Nrep = 3;

op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};

config = opt_config();
config.syn.elimination = false;


%form the system
bind =[ones(Nrep, 1); 2];
sys = opt_system(ops, [], [], bind);


%set the information structure

info_1 = tril(ones(Nrep), -1);

config.syn.D_mask = blkdiag(info_1, 1); 
man = opt_synthesis(sys, config);

Niter = 2;
order = {1, 1};
[sol_rep, history_rep] = man.alternate(Niter, order);


%% compare against an algorithm without repeated evaluation
sys_simple = opt_system(ops,[], []);
config_simple = opt_config;
config_simple.syn.elimination = false;
config_simple.syn.D_mask = [0, 0; 0, 1];

man_simple = opt_synthesis(sys_simple, config_simple);
[sol_simple, history_simple] = man_simple.alternate(Niter, order);


%% simulate and plot
d =  100;
Q1 = rand_quad(d, m, L); bstar1 = randn(d, 1)*100 - 20;
op1_sim = op_sim_quad(Q1, bstar1);

%define the L infinity ball (as a test)
BOX = 10;
op2_sim = op_sim_box(BOX);
ops_sim = {op1_sim, op2_sim};


sys_sim = sol_rep{end, end}.sys.export_sim(ops_sim);
sys_sim_simple = sol_simple{end, end}.sys.export_sim(ops_sim);
sim = alg_sim(sys_sim, d);
sim_simple = alg_sim(sys_sim_simple, d);
T = 30;

sim_out= sim.sim(T);
sim_out_simple= sim_simple.sim(T);


%% plot the result
figure(4)
clf
tiledlayout(1, 2)
nexttile
hold on
plot(sim_out_simple.k, sim_out_simple.res_w, 'LineWidth',3)
plot(sim_out.k, sim_out.res_w, '-.', 'LineWidth',3)
set(gca, 'yscale', 'log')
title('Optimality Error', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$||1^{\top} w||_2$', 'interpreter', 'latex', 'fontsize', 14)
xlabel('$k$', 'interpreter', 'latex', 'fontsize', 14)
legend({'$h=1$', '$h=3$'}, 'interpreter', 'latex', 'fontsize', 14)

nexttile 
hold on
plot(sim_out_simple.k, sim_out_simple.res_z, 'LineWidth',3)
plot(sim_out.k, sim_out.res_z, '-.', 'LineWidth',3)
set(gca, 'yscale', 'log')
title('Consensus Error', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$||z - z_{avg}||_2$', 'interpreter', 'latex', 'fontsize', 14)
xlabel('$k$', 'interpreter', 'latex', 'fontsize', 14)
legend({'$h=1$', '$h=3$'}, 'interpreter', 'latex', 'fontsize', 14)

% plt = alg_plotter(sim_out);
% plt_simple = alg_plotter(sim_out_simple);
% plt.plot_6f(1);
% plt_simple.plot_6f(2);
