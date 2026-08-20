rng(32, 'twister');

d = 200; %dimension of variable beta

%% describe the operators
%define the quadratic
m = 1; L = 5;
Q = rand_quad(d, m, L);
zstar = 100*(2*rand(d, 1) - 1);
op1 = op_sim_quad(Q, zstar);
h=3;

%define the L1 ball
tau = 100;
op2 = op_sim_l1_hard(tau);
ops = {op1, op2};

%% describe the network
%use a transfer function representation 
%to model the channel memory
z = tf('z', 1);
alpha = 0.4; %channel memory effect


P = [0, 0, 1, 0;
    0 , 0, 0, 1;
    1/(z-alpha)^h, 0, 0, 0;
    0, 1, 0, 0];

%partition the input and output channels
network = genplant(P);
network.nw = 2; network.nu = 2;
network.nz = 2; network.ny = 2;

%% describe the controller
b0 = -0.04;
b1 = -0.2;
b2 = -0.1;

ascale =  (1-alpha)^h;
K = ss([1], [ascale*b0, b0], [1; 1], [ascale*b2, 0; ...
    ascale*(b1+b2), b1],1);


%% form the system
sys = opt_system(ops, network, K);


%% simulate and plot
sim = alg_sim(sys, d);
T = 100;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot_4(1);




%% regulator equation tracking

sim_out_long= sim.sim(3*T); %get a more accurate solution to judge tracking

betastar = sim_out_long.z(end, :, end);
wend = sim_out_long.w(:, :, end);
dstar = [-betastar; wend(1:end-1, :, end)]; %the tracked reference


%compute the regulator equation solution
reg = regulator_lti(sys);
regcl = reg.check_regulator();
plt = plt.add_opt_sig(regcl, dstar);
plt.plot_4_err(2);    %plot the tracking error
plt.plot_4_sq_err(3); %plot the  squared norm of the tracking error