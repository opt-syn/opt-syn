%% describe the operators
m = 1; L = 5;

SML = false;

if SML
    op1 = op_sml(m, L);
    order = {[1, 1], [1,1]};
else
    op1 = op_gen();
    op1.monotone = m;
    op1.cocoercive = 1/L;
    order = {1, 1};
end

op2 = op_pcc();
ops = {op1, op2};

%% describe the network
%use a transfer function representation 
%to model the channel memory
z = tf('z', 1);
% alpha = 0.4; %channel memory effect
alpha = 1.2;
ascale = (2*(alpha+1));

P = [0, 0, z/(z+alpha), 0;
    0 , 0, 0, 1;
    z/(z+alpha), 0, 0, 0;
    0, 1, 0, 0];

%partition the input and output channels
network = genplant(P);
network.nw = 2; network.nu = 2;
network.nz = 2; network.ny = 2;

%add performance requirement (l2 stability)
network = network.add_oracle_shift(1);

%define the specifications
spec = spec_l2(1);
spec.target = true;
specs = {spec};

%% form the system
sys = opt_system(ops, network);

%% solve the problem
config =opt_config();
config.syn.elimination = true;
config.gen.same_rho = true;
man = opt_synthesis(sys, config);


%three rounds of alternation
Niter = 3; 
[sol_h, v_h] = man.alternate(Niter, order);
sol = sol_h{end, end};
% sol = man.bisect([], specs);
% sol.rho


%% simulate and plot

rng(32, 'twister');

d = 200; %dimension of variable beta

%form the operators
Q = rand_quad(d, m, L);
zstar = 100*(2*rand(d, 1) - 1);
op1_sim = op_sim_quad(Q, zstar);


%define the L1 ball
tau = 100;
op2_sim = op_sim_l1_hard(tau);
ops_sim = {op1_sim, op2_sim};

sys_sim = sol.sys.export_sim(ops_sim);
sim = alg_sim(sys_sim, d);
wp_eps = 0;
sim.sampler.wp = @(k, param) wp_eps*(2*rand(1, d)-1);

T = 61;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_6(1);


%% regulator equation tracking



sim_out_long= sim.sim(3*T); %get a more accurate solution to judge tracking

betastar = sim_out_long.z(end, :, end);
wend = sim_out_long.w(:, :, end);
dstar = [-betastar; wend(1:end-1, :, end)]; %the tracked reference


%find tracking error
plt = plt.add_opt_sig(sol.regcl, dstar);
plt.plot_4_err(2);    %plot the tracking error
plt.plot_4_sq_err(3); %plot the  squared norm of the tracking error
% 
% %% compare against the algorithm in the simulation example
% gamma = 0.4;
% lambda = 0.2;
% 
% K = ss([1], [-gamma*lambda, -gamma*lambda/(alpha+1)], ...
%     [1+alpha; 1], [0, 0; -gamma, -gamma/(alpha+1)],1);
% sys_given = sys_sim;
% sys_given.K = K;
% sys_given.op = sys.op;
% man_ana = opt_analysis(sys_given);
% order_ana = {[2, 2], [2, 2]};
% sol_ana = man_ana.bisect(order_ana);
% rho_ana = sol_ana.rho