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

%form the system
sys = opt_system(ops, network, []);

%synthesize an algorithm
config =opt_config();
config.syn.prox = [0, 0, 1];
man= opt_synthesis(sys, config);
sol= man.bisect();

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
T = 2000;

sim  = alg_sim(sys, d);
nx = size(sys.get_alg(1).A, 1);
% sim.sampler.x0 = 4*randn(nx, d);
sim_out = sim.sim(T);

plt = alg_plotter(sim_out);
plt.plot({'xc', 'w', 'xn', 'z'}, 1)
plt.plot({'f', 'res_w', 'res_z' }, 5)



%% regulator equations
sim_out_long = sim.sim(2*T);
zend = sim_out_long.z(:, :, end);
betastar = zend(1, :);
wend = sim_out_long.w(:, :, end);
dstar = [-betastar; wend(1:end-1, :)];
plt = plt.add_opt_sig(sol.regcl, dstar);

%% new
plt = plt.add_opt_sig(sol.regcl, dstar);
plt.plot_4_err(2);    %plot the tracking error
plt.plot_4_sq_err(3); %plot the  squared norm of the tracking error