%pose the operators
m = 1; L = 3;
ops = {op_sml(m, L)};

%form the network
n = struct('nw', 1, 'nz', 1, 'ny', 1, 'nu', 1);    
A = [1.2, 0; 0, -0.2];
B = [1, 0; 0, 1];
C = [0, 1; 1, 0];
D = [0, 0; 0, 2];

network = genplant(ss(A, B, C, D, 1), n);

sys = opt_system(ops, network, []);

config =opt_config();
config.syn.D_mask = 0;
config.gen.same_rho = true;
config.syn.reduced_order = false;
man= opt_synthesis(sys, config);
sol= man.bisect();

%% begin simulation (attempted)
d = 30;
Q = rand_quad(d, m, L);
bstar = randi(101, [d, 1]) - 50;
% bstar = zeros(d, 1);
op1 = op_sim_quad(Q, bstar);

ops_sim = {op1};

sys = sol.sys.export_sim(ops_sim);

T = 200;

sim  = alg_sim(sys, d);
nx = size(sys.get_alg(1).A, 1);
sim.sampler.x0 = 4*randn(nx, d);
sim_out = sim.sim(T);

plt = alg_plotter(sim_out);
plt.plot({'xc', 'w', 'res_w', 'xn', 'z', 'f' }, 1)

%% regulator equations
sim_out_long = sim.sim(20*T);
dstar = -sim_out_long.z(:, :, end);
plt = plt.add_opt_sig(sol.regcl, dstar);

%% new
plt = plt.add_opt_sig(sol.regcl, dstar);
plt.plot_4_err(2);    %plot the tracking error
plt.plot_4_sq_err(3); %plot the  squared norm of the tracking error