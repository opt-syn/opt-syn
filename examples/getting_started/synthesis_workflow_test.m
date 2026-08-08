%describe the operators 
m = 1; L = 50;
op1 = op_sml(m, L); %gradient of f
op2 = op_pcc();     %indicator function of L1 ball


%run the synthesis routine, use bisection to minimize the convergence rate
sys = opt_system({op1, op2});

man = opt_synthesis(sys); 
sol = man.bisect();

rho = sol.rho; % 0.8676


%% with time-delay dynamics
delay = [1, 0];
network = bridge_channel_delay(delay, delay);
sys = opt_system({op1, op2}, network);
man = opt_synthesis(sys);
sol_delay = man.bisect();  % 0.9860


%% simulate as a test
d = 500; %dimension of problem

%function f: random quadratic
M = rand_quad(d, m, L);
zstar = randi(101, [d, 1]) - 50;
op1_sim = op_sim_quad(M, zstar);


%convex set Z: l1 norm ball
tau = 100;
op2_sim = op_sim_l1_hard(tau);

%form simulation system
ops_sim = {op1_sim, op2_sim};
sys_sim = sol.sys.export_sim(ops_sim);

%execute algorithm
sim = alg_sim(sys_sim, d);
T = 150; ssim= sim.sim(T);

% plot the signals
plt = alg_plotter(ssim);
plt.plot({'x', 'w', 'res_w', 'f', 'z', 'res_z'}, 13)