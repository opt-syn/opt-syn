%% describe the operators
%define the quadratics
s = 4;
m = [1, -2, 0, 2]; L = [2, -1, inf, 10];

CHANNEL = true;

ops = cell(s, 1);
for i = 1:s
    ops{i} = op_sml(m(i), L(i));
end

pat_prox_seq = tril(ones(s));
pat_grad_seq = tril(ones(s), -1);
pat_grad_seq(3, 3) = 1;
pat_prox_par = eye(s);
pat_grad_par = zeros(s);
pat_grad_par(3, 3) = 1;

mask_pattern = {pat_prox_seq, pat_prox_par, pat_grad_seq, pat_grad_par};
% mask_pattern = {pat_prox_seq};

%% describe the network
network_1 = bridge_pass_through(2);
%use a transfer function representation 
%to model the channel memory
z = tf('z', 1);
% alpha = 0.4; %channel memory effect
alpha = 2;
% ascale = (2*(alpha+1));

P = [0,0,0,0,   1,0,0,0;
     0,0,0,0,   0,z/(z+alpha), 0, 0;
     0,0,0,0,   0,0, 1, 0;
     0,0,0,0,   0, -2/z,0,1;
     1,0,0,0,   0,0,0,0;
     0, z/(z+alpha), 0,0,   0,1,0,0;
     0,0,1,0,         0,0,0,0;
     0,-2/z,0,1,   0,0,0,0];

%partition the input and output channels
network = genplant(P);
network.nw = 4; network.nu = 4;
network.nz = 4; network.ny = 4;

%% form the system
sys_net = opt_system(ops, network);
sys = opt_system(ops, []);

%% solve the problem

Nmask = length(mask_pattern);
rho_pat = zeros(4, Nmask);
% for i = 1:Nmask
for i = 2:2
    mp= mask_pattern{i};
    config =opt_config();

    %configuration for the first Structure, Analysis, Synthesis paper
    config.gen.same_rho = true;
    config.syn.reduced_order=false;
    config.syn.D_mask = mask_pattern{i};
    config.tol.GX_max = 500;
    config.tol.GY_max = 500;
    config.bisect.val_range(2) = 10;
    man = opt_synthesis(sys, config);    

    % %three rounds of alternation
    Niter = 3; 
    order = 2;
    [sol_h, v_h] = man.alternate(Niter, order);
    sol = sol_h{end, end};
%     sol = man.bisect();
    rho_pat(1, i) = sol_h{1, 1}.rho;    
    rho_pat(3, i) = sol_h{1, end}.rho;    

    man_net = opt_synthesis(sys_net, config); 
%     sol_net = man_net.bisect();
    [sol_h, v_h] = man_net.alternate(Niter, order);
%     rho_pat(2, i) = sol_net.rho;   
    rho_pat(2, i) = sol_h{1, 1}.rho;    
    if isempty(sol_h{1, end})
        rho_pat(4, i) = rho_pat(2, i);    
    else
        rho_pat(4, i) = sol_h{1, end}.rho;    
    end
end

%% simulate and plot

rng(32, 'twister');

d = 50; %dimension of variable beta
ops_sim = cell(s, 1);
%form the operators
for i = 1:s
    if L(i) < inf
        Q = rand_quad(d, m(i), L(i));
        zstar = 100*(2*rand(d, 1) - 1);
        ops_sim{i}= op_sim_quad(Q, zstar);
    else

        %define the L1 ball
        tau = 100;       
        ops_sim{i} = op_sim_l1_hard(tau);
    end
end

sys_sim = sol_net.sys.export_sim(ops_sim);
sim = alg_sim(sys_sim, d);
T = 50;
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