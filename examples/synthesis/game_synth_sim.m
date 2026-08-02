rng(33, 'twister');

%define the Linear Quadratic Game
N_agent = 4; %number of agents
n = 5*ones(N_agent, 1); %strategy space of each agent
mu_boost = 5; %add strong monotonicity
[Q_list, b_list, c_list] = rand_LQ_game(n, mu_boost);
op1_sim = op_sim_LQ_game(Q_list, b_list, c_list, n);


%% synthesize an algorithm
mu = op1_sim.monotone; %1.4785
beta = op1_sim.coco; %0.1605

op1 = op_gen();
op1.monotone = mu;
op1.cocoercive = beta;

op2 = op_pcc();
ops = {op1, op2};

%form the network delay
d1 = [1, 0];
d2 = [1, 0];
network = bridge_channel_delay(d1, d2);

%form the system
sys = opt_system(ops, [], []);
sys_delay = opt_system(ops, network, []);

%set the information structure
config = opt_config();
config.syn.elimination = false;

%backward evaluation of pseudogradient
config.syn.D_mask = [1, 0; 1, 1]; 
man_bw = opt_synthesis(sys, config);

%forward evaluation of pseudogradient
config.syn.D_mask = [0, 0; 1, 1]; 
man_fw = opt_synthesis(sys, config);

%delay on pseudogradient, no controller restriction
config.syn.D_mask = [1, 1; 1, 1]; 
man_delay = opt_synthesis(sys_delay, config);

sol_bw = man_bw.bisect();       %rho = 0.7163
sol_fw = man_fw.bisect();       %rho = 0.8734
sol_delay = man_delay.bisect(); %rho = 0.9430

%% simulate and plot

%define the L infinity ball
%individual constraints for agents
BOX = 10;
op2_sim = op_sim_box(BOX);
ops_sim = {op1_sim, op2_sim};

sys_sim_fw = sol_fw.sys.export_sim(ops_sim);
sys_sim_bw = sol_bw.sys.export_sim(ops_sim);
sys_sim_delay = sol_delay.sys.export_sim(ops_sim);


sys_sim_list = {sys_sim_fw, sys_sim_bw, sys_sim_delay};

for i =1:3
    d = sum(n);
    sim = alg_sim(sys_sim_list{i}, d);
    T = 100;
    
    sim_out= sim.sim(T);
    plt = alg_plotter(sim_out);
    plt.plot({'x', 'w', 'res_w', 'payoff', 'z', 'res_z'}, i);
end