%channel memory

m =1;
L = 5;
%get the operator

op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};
m = 1; L = 10;



gamma = 0.4;
lambda = 0.2;

%use a transfer function representation 
%to model the channel memory
z = tf('z', 1);
alpha = 0.4;
ascale = (2*(alpha+1));

P = [0, 0, z/(z+alpha), 0;
    0 , 0, 0, 1;
    z/(z+alpha), 0, 0, 0;
    0, 1, 0, 0];

network = genplant(P);
network.nw = 2; network.nu = 2;
network.nz = 2; network.ny = 2;

K = ss([1], [-gamma*lambda, -gamma*lambda/(alpha+1)], ...
    [1+alpha; 1], [0, 0; -gamma, -gamma/(alpha+1)],1);

sys = opt_system(ops, network, K);

reg = regulator_lti(sys);
regcl = reg.check_regulator();



man = opt_analysis(sys);
order = {1, 1};
sol = man.bisect(order);
sol.rho
