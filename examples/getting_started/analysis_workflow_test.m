% The set of considered oracles
m = 1;
L = 10;
op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};

% We solve the optimization problem with a Douglas-Rachford algorithm 
%Douglas-Rachford algorithm
% gamma = 0.4;
% lambda = 0.25;
% K = ss(1, [-gamma*lambda, -gamma*lambda], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);

%PGD
gamma = 2/(L + m);
K = ss(1, [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);


%form the system and the analysis manager
sys = opt_system(ops, [], K);
man = opt_analysis(sys);

% Define the order 
% order = {[1, 1], [1, 1]};
order = {1, 1};
% order = {[2, 2], [2, 2]};
sol_best = man.bisect(order);

% %sweep
% nL = 300;
% Lsweep = linspace(1, 5000, nL);
% rholist = zeros(nL, 1);
% for i = 1:length(Lsweep)
%     man.sys.op{1}.L = Lsweep(i);
%     sol_best = man.bisect(order);
%     rholist(i) = sol_best.rho;
% end

%% delay
% We introduce a time delay of 2 steps before and after the first oracle evaluation
DELAY = [2, 0];
network_delay = bridge_channel_delay(DELAY, DELAY);

% order_delay = {[2, 2], [2, 2]};
order_delay = {1, 1};
sys_delay = opt_system(ops, network_delay, K);
man_delay = opt_analysis(sys_delay);
[sol_best_delay, v_range_delay] = man_delay.bisect(order_delay);