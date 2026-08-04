%coordinate descent scheme
rng(50, 'twister');

d = 90; %number of dimensions
% c=6; %blocksize 
c = 3;

%symmetry generator/permutation matrix
M = circshift(eye(c), 1); 

%define the objective function
m = 1; L = 5; 

op1 = op_sml(1, 5);

op1.c = c; %enforce coordinate dimension

%define the network
%implements coordinate descent

network = coordinate_descent_primitive(c);
%form the system
sys = opt_system_periodic_orbit(op1, network, [], M);


config = opt_config();
config.syn.D_mask = 0;

man = opt_synthesis(sys, config);

% %access to all parts of gradient
% sol = man.bisect();


%only access to a part of the gradient

%lowered from default
config.tol.spread = 1e-4;
config.tol.input_diss = 1e-5;
config.tol.M = 1e-9;

%raised from default
config.tol.GX_max = 300;   
config.tol.GY_max = 300;


network_lim = coordinate_descent_limited(c);
sys_lim = opt_system_periodic_orbit(op1, network_lim, [], M);
man_lim = opt_synthesis(sys_lim, config);
[sol_alt, v_history] = man.alternate(1, {2}, []);
% sol = man_lim.bisect();
% sol = sol_alt(1);

%% simulate coordinate descent
%define the objective function
%sum of quadratic and log-sum-exp
L2 = 3;
Q = rand_quad(d, m, L2);
bstar = 100*(2*rand(d, 1) - 1);

fw_func = @(k, z, param) Q*z + bstar + (L - L2)* (exp(z) - exp(-z)) / sum(exp(z) + exp(-z));
bw_func = @(k, z, D, param) []; %unused
f_func =  @(k, z, param) 0.5*z'*Q*z + bstar' * z + (L - L2)*log(sum(exp(z) + exp(-z)));
op1_sim = op_sim(fw_func, bw_func, f_func);

ops_sim= {op1_sim};
sys_sim = sol.sys.export_sim(ops_sim);
sims = alg_sim(sys_sim, d);



sys_lim_sim = sol_lim.export_sim(ops_sim);
sims_lim = alg_sim(sys_lim_sim, d);

T = 50;
sim_out= sims.sim(T);
sim_out_long= sims.sim(50*T);

T_lim = 1000;
sim_lim_out= sims_lim.sim(T_lim);


%% plot
plt = alg_plotter(sim_out);
plt.plot({'xn', 'w', 'res_w', ...
    'xc', 'z', 'coord', }, 10);