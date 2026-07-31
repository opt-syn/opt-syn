rng(33, 'twister');


DR = true;
if DR
    %Douglas Rachford Algorithm
    gamma = 1;
    lambda = 1;

    K = ss([1], [-lambda*gamma, -lambda*gamma], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);

else

    %projected pseudogradient descent
    gamma = 0.9;    %stepsizes
    sK0 = ss([1], [-gamma, -gamma], [1; 1], [0,  0; -gamma, -gamma],1);
end

%define the game
N_agent = 4;
n = 5*ones(N_agent, 1); %strategy space of each agent
mu_boost = 5;
[Q_list, b_list, c_list] = rand_LQ_game(n, mu_boost);
op1 = op_sim_LQ_game(Q_list, b_list, c_list, n);


%define the L infinity ball
BOX = 10;
op2 = op_sim_box(BOX);
ops = {op1, op2};

%douglas-rachford
gamma = 0.3;    %stepsizes
% K = ss(1, [-gamma, -gamma], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);
K = ss(1, [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);

%form the system
sys = opt_system(ops, [], K);

%% simulate and plot
d = sum(n);
sim = alg_sim(sys, d);
T = 30;
% sim.sampler.x0 = 200*randn(1, d);
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot({'x', 'w', 'res_w', 'payoff', 'z', 'res_z'});


%% regulator equation tracking
sim_out_long= sim.sim(5*T);
zend = sim_out_long.z(:, :, end);

betastar = sim_out_long.z(end, :, end);
wend = sim_out_long.w(:, :, end);
dstar = [-betastar; wend(1:end-1, :, end)];

reg = regulator_lti(sys);
regcl = reg.check_regulator();
plt = plt.add_opt_sig(regcl, dstar);
plt.plot_3_err(2);