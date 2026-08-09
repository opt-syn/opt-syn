%pose the operators
m = 1; L = 3;
ops = {op_sml(m, L), op_sml(1, 1), op_pcc()};

%form the network
n = struct('nw', 1, 'nz', 1, 'ny', 1, 'nu', 1);    
A = [1.2, 0; 0, -0.2];
B = [1, 0; 0, 1];
C = [0, 1; 1, 0];
D = [0, 0; 0, 2];
P = ss(A, B, C, D, 1);
network1 = genplant(P, n);
network2 = bridge_pass_through(2);
network = blkdiag(network1, network2);

%add robustness (l2 gain)
%modify the network
[network, iwp] = network.add_oracle_input(1, []);
[network, izp] = network.perf_output_con();
%impose specification
GAIN = 10;
spec = spec_e2e(GAIN, iwp, izp);
spec.target = true;
specs = {spec_stability(1), spec};

%form the system
sys = opt_system(ops, network, []);


%synthesize an algorithm
config =opt_config();
config.syn.prox = [0, 0, 1];
config.bisect.bisect_rho = true;
% config.syn.reduced_order = false;
config.gen.same_rho = true;
man= opt_synthesis(sys, config);
sol= man.bisect([], specs);

%% begin simulation (attempted)
d = 100;
Q = rand_quad(d, m, L);
bstar = randi(101, [d, 1]) - 50;
% bstar = zeros(d, 1);
op1 = op_sim_quad(Q, bstar);


bstar2 = randi(101, [d, 1]) + 50;
op2 = op_sim_quad(eye(d), bstar2);

tau = 50;
op3 = op_sim_l1_hard(tau);

ops_sim = {op1, op2, op3};

sys = sol.sys.export_sim(ops_sim);

% T = 301;
T = 600;
% T = 2000;

sim  = alg_sim(sys, d);
eps_w = 0;
sim.sampler.wp = @(k, param) (2*rand(1, d)-1)*eps_w;
nx = size(sys.get_alg(1).A, 1);
% sim.sampler.x0 = 4*randn(nx, d);
sim_out = sim.sim(T);


%% regulator equations
sim_long = sim;
sim_long.sampler.wp = @(k, param) (2*rand(1, d)-1)*0;
sim_out_long = sim_long.sim(2*T);

zend = sim_out_long.z(:, :, end);
betastar = zend(1, :);
wend = sim_out_long.w(:, :, end);
dstar = [-betastar; wend(1:end-1, :)];

%% new

plt = alg_plotter(sim_out);
plt = plt.add_opt_sig(sol.regcl, dstar);

plt.plot({'xc', 'w', 'xn', 'z'}, 1)
plt.plot({'sq_xerr', 'res_w', 'res_z' }, 5)

plt.plot_4_err(2);    %plot the tracking error
plt.plot_4_sq_err(3); %plot the  squared norm of the tracking error