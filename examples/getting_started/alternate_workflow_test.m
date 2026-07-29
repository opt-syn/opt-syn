% The set of considered oracles
m = 1;
L = 10;
op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};
sys = opt_system(ops);
man = opt_synthesis(sys);

%% bisection
sol_bisect= man.bisect();

%% alternation
Niter = 3;
[sol_alt, v_r] = man.alternate(Niter,  {1, 1});

%% simulate as a test
sol = sol_alt{1, end};

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