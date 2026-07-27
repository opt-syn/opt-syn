% The set of considered oracles
m = 1;
L = 10;
% L = 3;
op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};

%form the system and the analysis manager
sys = opt_system(ops, [], []);
config = opt_config();
config.syn.reduced_order = true;
config.syn.elimination = true;
% config.syn.reduced_order = ;
% config.syn.elimination = false;
man = opt_synthesis(sys, config);
sol_best = man.bisect();

% spec = spec_stability(0.5);
% sol_single = man.solve_single([], spec)
% 
% %% gradient
% config_grad = opt_config();
% config_grad.syn.D_mask = [0, 0; 1, 1];
% man_grad = opt_synthesis(sys, config_grad);
% sol_grad = man_grad.bisect();
% 
% %% alternation
% man_alt = man;
% Niter =3;
% sol_alt = man_alt.alternate(3, {1, 1});

%% simulate as a test
sol = sol_best;

d = 50;
BOX = 30;
M = rand_quad(d, m, L);
zstar = randi(101, [d, 1]) - 50;

op1_sim = op_sim_quad(M, zstar);
op2_sim = op_sim_box(BOX);

ops_sim = {op1_sim, op2_sim};

sys_sim = sol.sys.export_sim(ops_sim);

sim = alg_sim(sys_sim, d);
T = 100;
ssim= sim.sim(T);

% plot the signal
plt = alg_plotter(ssim);
plt.plot({'x', 'w', 'res_w', 'f', 'z', 'res_z'}, 13)