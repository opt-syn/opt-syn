rng(32, 'twister');

%no network dynamics
%Douglas Rachford Algorithm

d = 300; %dimension of variable beta


%douglas-rachford
gamma = 0.4;    %stepsizes
lambda = 0.25;
sK0 = ss([1], [-gamma*lambda, -gamma*lambda], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);


%define the quadratic
m = 1; L = 1000;
Q = rand_quad(d, m, L);
zstar = 100*(2*rand(d, 1) - 1);
op1 = op_sim_quad(Q, zstar);


%define the L infinity ball
BOX = 10;
op2 = op_sim_box(BOX);
ops = {op1, op2};

%douglas-rachford
gamma = 0.4;    %stepsizes
lambda = 0.25;
sK0 = ss(1, [-gamma*lambda, -gamma*lambda], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);

%form the system
sys = opt_system(ops, [], K);

%% simulate and plot
sim = alg_sim(sys, d);
T = 20;
sim.sampler.x0 = 200*randn(1, d);
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
fig1 = plt.plot_4(1);


%% regulator equation tracking
zend = sim_out.z(:, :, end);

betastar = sim_out.z(end, :, end);
wend = sim_out.w(:, :, end);
dstar = [-betastar; wend(1:end-1, :, end)];
reg = regulator_lti(sys);
regcl = reg.check_regulator();
plt = plt.add_opt_sig(regcl, dstar);
fig2 = plt.plot_3_sq_err(2);
