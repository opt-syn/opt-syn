c=6; %blocksize 

%symmetry generator/permutation matrix
M = circshift(eye(c), 1); 

%define the objective function
m = 1; L = 2; 

op1 = op_sml(m, L);

op1.c = c; %enforce coordinate dimension

%performs gradient operations, not prox
config = opt_config();
config.syn.D_mask = 0;

%% Full Knowledge
network = coordinate_descent_system(c);

%form the system
sys = opt_system_periodic_orbit(op1, network, [], M);

man = opt_synthesis(sys, config);

% %access to all parts of gradient
sol = man.bisect();

%% Partial Knowledge
network_lim = coordinate_descent_limited(c);
sys_lim = opt_system_periodic_orbit(op1, network_lim, [], M);
man_lim = opt_synthesis(sys_lim, config);
sol_lim = man_lim.bisect();


%compare the solutions
[sol.rho, sol_lim.rho] 

%c=1: 0.3557    0.3557
%c=2: 0.5829    0.7487
%c=3: 0.6958    0.8750 
%c=4: 0.7640    0.9304
%c=5: 0.8108    0.9595
%c=6: 0.8441    0.9767

%% simulate coordinate descent
%define the objective function
%sum of quadratic and log-sum-exp
d = 120; %number of dimensions
L2 = 1.5;
rng(50, 'twister');


Q = rand_quad(d, m, L2);
bstar = 100*(2*rand(d, 1) - 1);

fw_func = @(k, z, param) Q*z + bstar + (L - L2)* (exp(z)) / sum(exp(z));
bw_func = @(k, z, D, param) []; %unused
f_func =  @(k, z, param) 0.5*z'*Q*z + bstar' * z + (L - L2)*log(sum(exp(z)));
op1_sim = op_sim(fw_func, bw_func, f_func);

ops_sim= {op1_sim};

% execute the algorithms
sys_sim = sol.sys.export_sim(ops_sim);
sys_lim_sim = sol_lim.sys.export_sim(ops_sim);

sims = alg_sim(sys_sim, d);
sims_lim = alg_sim(sys_lim_sim, d);

T = 100;
sim_out= sims.sim(T);
sim_lim_out= sims_lim.sim(T);


%% plot
sigs = {'xn', 'w', 'res_w', 'xc', 'z', 'coord'};
plt = alg_plotter(sim_out);
plt.plot(sigs, 10);

plt_lim = alg_plotter(sim_lim_out);
plt_lim.plot(sigs, 11);


%% plot c
rho = [0.3557    0.3557;
    0.5829    0.7487;
    0.6958    0.8750;
     0.7640    0.9304;
    0.8108    0.9595;
     0.8441    0.9767;
       0.8686    0.9877;
         0.8874    0.9949];
figure(40)
clf
Nc = size(rho, 1);
scatter(1:Nc, rho, 30, '-o', 'Linewidth', 2, 'filled');
xlabel('$c$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\rho$', 'interpreter', 'latex', 'fontsize', 16)
legend({'Full Knowledge', 'Partial Knowledge'}, 'location','southeast', 'interpreter', ...
    'latex', 'fontsize', 16)
ylim([0, 1])