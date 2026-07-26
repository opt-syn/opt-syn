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
order = {1, 1};

%get the convergence rate
rho = 0.95;
spec = spec_stability(rho);

sol_best = man.solve_single(order, spec)
