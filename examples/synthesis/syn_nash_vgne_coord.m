%synthesis of fixed point algorithms 
%set-valued maps (nash equilibria)

%% generate parameters of the game
rng(40, 'twister')

c=4;

mu = 1.4785;
beta = 0.1605;

BOX = 5;

op1 = op_gen();
op1.monotone = mu;
op1.cocoercive = beta;

op2 = op_pcc();
ops = {op1, op2};

M =  circshift(eye(c), -1);

network = blkdiag(coordinate_descent_system(c), bridge_pass_through(c));
% network = blkdiag(coordinate_descent_limited(c), bridge_pass_through(c));

sys_simul = opt_system(ops);
sys_coord = opt_system_periodic_orbit(ops,  network, [], M);

config = opt_config();
config.syn.prox= [0, 1];

man_simul = opt_synthesis(sys_simul, config);
man_coord = opt_synthesis(sys_coord, config);


sol_simul = man_simul.bisect();
sol_coord = man_coord.bisect();

%% plot

rng(33, 'twister');

%define the  game
%each agent minimizes a quadratic cost
N_agent = 4; %number of agents
n = 5*ones(N_agent, 1); %strategy space of each agent
mu_boost = 5; %add strong monotonicity
[Q_list, b_list, c_list] = rand_LQ_game(n, mu_boost);
op1_sim = op_sim_LQ_game(Q_list, b_list, c_list, n);


sys_sim = sol_simul.sys.export_sim({op1_sim});
sys_sim_coord = sol_coord.sys.export_sim({op1_sim});

% simulate and plot
d = sum(n);
sim = alg_sim(sys_sim, d);
sim_coord = alg_sim(sys_sim_coord, d);
T = 61;
% sim.sampler.x0 = 200*randn(1, d);
sim_out= sim.sim(T);
sim_out_coord= sim_coord.sim(T);

%% plot the results
plt = alg_plotter(sim_out);
% plt.plot({'x', 'w', 'res_w', 'payoff', 'z', 'coord'}, 1);
plt.plot({'x', 'res_w', 'payoff', 'z'}, 1);


plt_coord = alg_plotter(sim_out_coord);
plt_coord.plot({'x', 'res_w', 'payoff', 'z'}, 2);
