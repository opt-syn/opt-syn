rng(33, 'twister');

%define the  game
%each agent minimizes a quadratic cost
N_agent = 4; %number of agents
n = 5*ones(N_agent, 1); %strategy space of each agent
mu_boost = 5; %add strong monotonicity
[Q_list, b_list, c_list] = rand_LQ_game(n, mu_boost);
op1 = op_sim_LQ_game(Q_list, b_list, c_list, n);


%define the L infinity ball
%individual constraints for agents
BOX = 10;
op2 = op_sim_box(BOX);
% op2 = op_sim_l1_hard(BOX);
ops = {op1, op2};

%douglas-rachford
gamma = 1;
lambda = 1;

K = ss([1], [-lambda*gamma, -lambda*gamma], ...
    [1; 1], [-gamma, 0; -2*gamma, -gamma],1);

%form the system
sys = opt_system(ops, [], K);

%% simulate and plot
d = sum(n);
sim = alg_sim(sys, d);
T = 100;
% sim.sampler.x0 = 200*randn(1, d);
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot({'x', 'w', 'res_w', 'payoff', 'z', 'res_z'}, 1);