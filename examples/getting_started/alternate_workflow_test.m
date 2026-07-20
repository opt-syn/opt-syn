% The set of considered oracles
m = 1;
L = 10;
op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};

%% alternation

% config_grad = opt_config();
% % config_grad.syn.D_mask = [0, 0; 1, 1];
% man_alt = man_grad;
man_alt = man;
b_opts = bisect_opts;
b_opts.Niter = 3;
[sol_alt, v_r] = man_alt.alternate([], {1, 1}, [], b_opts);

%% simulate as a test
sol = sol_alt{1, end};

d = 50;
BOX = 30;
M = rand_quad(d, m, L);
zstar = randi(101, [d, 1]) - 50;

op1_sim = op_sim(@(k, z,  param) M*(z-zstar),...
    @(k, D, z, param) (M + kron(inv(D), eye(d))) \ (M*zstar + kron(D, eye(d)) \ z),...
    @(k, z, param) 0.5*(z-zstar)'*M*(z-zstar));

op2_sim = op_sim(@(k, z, param) zeros(size(z)),...
    @(k, D, z, param) clip(z, -BOX, BOX),...
    @(k, z, param) []);

ops_sim = {op1_sim, op2_sim};

sys_sim = opt_system(ops_sim, [], sol.K);

sim = alg_sim(sys_sim, d);
T = 100;
ssim= sim.sim(T);

%% plot the signal
plt = alg_plotter(ssim);
plt.plot({'x', 'w', 'res_w', 'f', 'z', 'res_z'}, 13)