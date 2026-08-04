%describe the operators 
m = 1; L = 50;
op1 = op_sml(m, L); %gradient of f
op2 = op_pcc();     %indicator function of L1 ball


%run the synthesis routine, use bisection to minimize the convergence rate
sys = opt_system({op1, op2});

man = opt_synthesis(sys); 
sol = man.bisect();

rho = sol.rho; % 0.9050


%% with time-delay dynamics
delay = [1, 0];
network = bridge_channel_delay(delay, delay);
sys = opt_system({op1, op2}, network);
man = opt_synthesis(sys);
sol = man.bisect();


%% simulate as a test
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