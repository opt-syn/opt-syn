%simulation of a multi-step algorithm
rng(40, 'twister');
m = [1; 0];
L = [3; 6];

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
d = 200; %dimension of problem
BOX = 30; %L infinity norm constraint

op1 = op_sml(m(1), L(1));
op2 = op_sml(m(2), L(2));
op3 = op_pcc();
ops = {op1, op2, op3};

%put the system together
bind = [1, 2, 3, 2];
sys = opt_system(ops, [], K, bind);

reg = regulator_lti(sys);
regcl = reg.check_regulator();


man = opt_analysis(sys);
order = {2, 2, 2};
sol = man.bisect(order)