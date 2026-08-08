% The set of considered oracles
m = 1;
L = 10;
op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};

%PGD
gamma = 2/(m+L);
K = ss(1, [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);


%form the system and the analysis manager
sys = opt_system(ops, [], K);
man = opt_analysis(sys);

% Define the order 
order = {1, 1};

%get the convergence rate
% rho = 0.95;
% spec = spec_stability(rho);

sol_best = man.bisect(order) %0.8182

%% with a time delay
delay = [1, 0];
network = bridge_channel_delay(delay, delay);
sys_delay = opt_system({op1, op2}, network, K);
man_delay = opt_analysis(sys_delay);
sol_delay = man_delay.bisect(order);  % 1.3744
