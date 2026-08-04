%coordinate descent scheme
rng(50, 'twister');

d = 30; %number of dimensions
% c=6; %blocksize 
% c = 3;
c = 5;

%symmetry generator/permutation matrix
M = circshift(eye(c), 1); 

%gradient descent rule
gamma = 0.2;

alg = ss(eye(c), blkdiag(-gamma , zeros(c-1)), ...
    eye(c), zeros(c), 1);

% K = ss(eye(c), blkdiag(-gamma , zeros(c-1)), ...
%     blkdiag(1, zeros(c-1)), zeros(c), 1);

% K = ss(eye(c), -gamma*eye(c), ...
%     eye(c), zeros(c), 1);


% K = ss(eye(c), blkdiag(-gamma , zeros(c-1)), ...
%     blkdiag(1, zeros(c-1)), zeros(c), 1);

%define the objective function
%sum of quadratic and log-sum-exp
m = 1; L = 5; L2 = 3;
Q = rand_quad(d, m, L2);
bstar = 100*(2*rand(d, 1) - 1);

fw_func = @(k, z, param) Q*z + bstar + (L - L2)* (exp(z) - exp(-z)) / sum(exp(z) + exp(-z));
bw_func = @(k, z, D, param) []; %unused
f_func =  @(k, z, param) 0.5*z'*Q*z + bstar' * z + (L - L2)*log(sum(exp(z) + exp(-z)));
op1 = op_sim(fw_func, bw_func, f_func);

op1.c = c; %enforce coordinate dimension

%define the network
%implements coordinate descent
% network = coordinate_descent_limited(c);

%form the system
% sys = opt_system_periodic_orbit(op1, network, K, M);
sys = opt_system_periodic_orbit(op1, [], alg, M);
reg = regulator_periodic_orbit(sys);
rcl = reg.check_regulator();

%% simulate coordinate descent
sim = alg_sim(sys, d);
T = 51;
sim_out= sim.sim(T);

% T_long = 801;
% sim_out_long= sim.sim(T_long);


%% plot
plt = alg_plotter(sim_out);
plt.plot({'f', 'w', 'res_w', ...
    'x', 'z', 'coord', }, 10);

% plt_long = alg_plotter(sim_out_long);
% plt_long.plot({'xn', 'w', 'res_w', ...
%     'xc', 'z', 'coord', }, 2);