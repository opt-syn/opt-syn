%coordinate descent scheme
rng(50, 'twister');

d = 90; %number of dimensions
% c=6; %blocksize 
c = 3;

%symmetry generator/permutation matrix
M = circshift(eye(c), 1); 

%define the objective function
%sum of quadratic and log-sum-exp
m = 1; L = 5; L2 = 3;
% Q = rand_quad(d, m, L2);
% bstar = 100*(2*rand(d, 1) - 1);
% 
% fw_func = @(k, z, param) Q*z + bstar + (L - L2)* (exp(z) - exp(-z)) / sum(exp(z) + exp(-z));
% bw_func = @(k, z, D, param) []; %unused
% f_func =  @(k, z, param) 0.5*z'*Q*z + bstar' * z + (L - L2)*log(sum(exp(z) + exp(-z)));
% op1 = op_sim(fw_func, bw_func, f_func);
op1 = op_sml(1, 5);

op1.c = c; %enforce coordinate dimension

%define the network
%implements coordinate descent
% network = coordinate_descent_limited(c);
Pprim = coordinate_descent_primitives(c);
network  = Pprim{1};
%form the system
sys = opt_system_periodic_orbit(op1, network, [], M);


config = opt_config();
config.syn.D_mask = 0;
man = opt_synthesis(sys, config);

sol = man.bisect();
% solh = man.alternate(3, 1, []);
% sol = solh{1, end};




%% simulate coordinate descent
%define the objective function
%sum of quadratic and log-sum-exp
m = 1; L = 5; L2 = 3;
Q = rand_quad(d, m, L2);
bstar = 100*(2*rand(d, 1) - 1);

fw_func = @(k, z, param) Q*z + bstar + (L - L2)* (exp(z) - exp(-z)) / sum(exp(z) + exp(-z));
bw_func = @(k, z, D, param) []; %unused
f_func =  @(k, z, param) 0.5*z'*Q*z + bstar' * z + (L - L2)*log(sum(exp(z) + exp(-z)));
op1_sim = op_sim(fw_func, bw_func, f_func);

ops_sim= {op1_sim};
sys_sim = sol.sys.export_sim(ops_sim);
sims = alg_sim(sys_sim, d);
T = 40;
sim_out= sims.sim(T);

T_long = 801;
sim_out_long= sims.sim(T_long);


%% plot
plt = alg_plotter(sim_out);
plt.plot({'xn', 'w', 'res_w', ...
    'xc', 'z', 'coord', }, 1);

plt_long = alg_plotter(sim_out_long);
plt_long.plot({'xn', 'w', 'res_w', ...
    'xc', 'z', 'coord', }, 2);