%simulation of a multi-step algorithm
rng(40, 'twister');
d = 200; %dimension of problem



%% create a multi-step controller that satisfies the regulator equation
gamma = 0.4;    %stepsizes
lambda = 0.25;


A = 1;
B = -gamma*lambda * [2, 1, 2, 1];
C = ones(4, 1);
D = -lambda * [1, 0, 0, 0;
               1, 0, 0, 0;
               2, 1, 1, 0;
               2, 1, 1, 0];


K = ss(A, B, C, D, 1);


%% form the operators


%randomly generate the quadratics
m1 = 1; m2 = 0;
L1 = 3; L2 = 6;
Q1 = rand_quad(d, m1, L1); bstar1 = randn(d, 1)*100 - 20;
Q2 = rand_quad(d, m2, L2); bstar2 = randn(d, 1)*100 + 20;

BOX = 30; %L infinity norm constraint

%create the operators
op1 = op_sim_quad(Q1, bstar1);
op2 = op_sim_quad(Q2, bstar2);
op3 = op_sim_box(BOX);
ops = {op1, op2, op3};

%put the system together
bind = [1, 2, 3, 2]; %assign the ordering of operators
sys = opt_system(ops, [], K, bind);

%% simulate and plot
sim = alg_sim(sys, d);
% T = 200;
T = 100;
sim.sampler.x0 = 200*randn(1, d);
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);
fig1 = plt.plot_4(1);

zend = sim_out.z(:, :, end);