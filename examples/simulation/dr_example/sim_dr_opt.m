rng(32, 'twister');

%no network dynamics
%Douglas Rachford Algorithm

d = 100; %dimension of variable beta




%define the quadratic
m = 1; L = 10;
Q = rand_quad(d, m, L);
zstar = 100*(2*rand(d, 1) - 1);
op1 = op_sim_quad(Q, zstar);


%douglas-rachford
gamma = 0.4;
lambda = 1;
K = ss([1], [-gamma*lambda, -gamma*lambda], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);


%define the L infinity ball
BOX = 10;
op2 = op_sim_box(BOX);
ops = {op1, op2};


%form the system
sys = opt_system(ops, [], K);

%% simulate and plot
sim = alg_sim(sys, d);
T = 100;
sim.sampler.x0 = 200*randn(1, d);
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
fig1 = plt.plot_6f(1);


%% regulator equation tracking

sim_out_long= sim.sim(10*T);
betastar = sim_out_long.z(end, :, end);
wend = sim_out_long.w(:, :, end);
dstar = [-betastar; wend(1:end-1, :, end)];
reg = regulator_lti(sys);
regcl = reg.check_regulator();
plt = plt.add_opt_sig(regcl, dstar);
fig2 = plt.plot_3_sq_err(2);
