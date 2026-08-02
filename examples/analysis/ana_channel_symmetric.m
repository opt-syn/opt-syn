%% describe the operators
%define the quadratic
m = 1; L = 5;

op1 = op_sml(m, L);
op2 = op_pcc();
ops ={op1, op2};

%% describe the network
%use a transfer function representation 
%to model the channel memory
z = tf('z', 1);
alpha = 0.4; %channel memory effect
ascale = (2*(alpha+1));

P = [0, 0, z/(z+alpha), 0;
    0 , 0, 0, 1;
    z/(z+alpha), 0, 0, 0;
    0, 1, 0, 0];

%partition the input and output channels
network = genplant(P);
network.nw = 2; network.nu = 2;
network.nz = 2; network.ny = 2;

%% describe the controller
gamma = 0.4;
lambda = 0.2;
K = ss([1], [-gamma*lambda, -gamma*lambda/(alpha+1)], ...
    [1+alpha; 1], [0, 0; -gamma, -gamma/(alpha+1)],1);


%% form the system
sys = opt_system(ops, network, K);

%% describe the operators
man_ana = opt_analysis(sys);
% order_ana = {[2, 2], [2, 2]};
order_ana = {1, 1};
sol_ana = man_ana.bisect(order_ana);
sol_ana.rho