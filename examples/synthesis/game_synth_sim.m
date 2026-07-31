rng(33, 'twister');


%form the network
d1 = [1, 0];
d2 = [1, 0];
network = bridge_channel_delay(d1, d2);


%define the noncooperative game
%each agent minimizes a quadratic cost
N_agent = 4; %number of agents
n = 5*ones(N_agent, 1); %strategy space of each agent
mu_boost = 5; %add strong monotonicity
[Q_list, b_list, c_list] = rand_LQ_game(n, mu_boost);
op1_sim = op_sim_LQ_game(Q_list, b_list, c_list, n);



%define the L infinity ball
%individual constraints for agents
BOX = 10;
op2_sim = op_sim_box(BOX);
ops_sim = {op1_sim, op2_sim};


%% synthesize an algorithm

mu = op1_sim.monotone;
beta = op1_sim.coco;

op1 = op_gen();
op1.monotone = mu;
op1.cocoercive = beta;

op2 = op_pcc();
ops = {op1, op2};


%form the system
sys = opt_system(ops, network, []);

config = opt_config();
config.syn.D_mask = [0, 0; 1, 1]; 

man = opt_synthesis(sys, config);
[sol_best, v_range] = man.bisect();
sol = sol_best;



%% simulate and plot



sys_sim = sol.sys.export_sim(ops_sim);

d = sum(n);
sim = alg_sim(sys_sim, d);
T = 100;
% sim.sampler.x0 = 200*randn(1, d);
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot({'x', 'w', 'res_w', 'payoff', 'z', 'res_z'}, 1);