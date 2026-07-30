rng(33, 'twister');

%Douglas Rachford Algorithm

%projected pseudogradient descent
gamma = 0.4;    %stepsizes
sK0 = ss([1], [-gamma, -gamma], [1; 1], [0,  0; -gamma, -gamma],1);


%define the game
N_agent = 4;
n = 5*ones(N_agent, 1); %strategy space of each agent
mu_boost = 27;
[Q_list, b_list, c_list] = rand_LQ_game(n, mu_boost);
op1 = op_sim_LQ_game(Q_list, b_list, c_list, n);


%define the L infinity ball
BOX = 5;
op2 = op_sim_box(BOX);
ops = {op1, op2};

%douglas-rachford
gamma = 0.2;    %stepsizes
K = ss(1, [-gamma*lambda, -gamma*lambda], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);

%form the system
sys = opt_system(ops, [], K);

%% simulate and plot
d = sum(n);
sim = alg_sim(sys, d);
T = 200;
sim.sampler.x0 = 200*randn(1, d);
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_4(1);


% %% regulator equation tracking
% zend = sim_out.z(:, :, end);
% 
% betastar = sim_out.z(end, :, end);
% wend = sim_out.w(:, :, end);
% dstar = [-betastar; wend(1:end-1, :, end)];
% reg = regulator_lti(sys);
% regcl = reg.check_regulator();
% plt = plt.add_opt_sig(regcl, dstar);
% plt.plot_3_err(2);