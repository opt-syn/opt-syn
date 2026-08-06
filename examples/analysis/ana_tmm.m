
%proximal triple-momentum method from https://arxiv.org/abs/2605.22929
%
%has convergence rate (1 - sqrt(m/L))

%define the operators
m = 1;
L = 10;

%different operator classes for op1
op1= op_sml(m, L);
op2 = op_sml(0, inf);
ops = {op1, op2};

%define the controller

A = [(1.0 - sqrt(m / L)) / (1.0 + sqrt(m / L)),  2.0 * sqrt(m / L) / (1.0 + sqrt(m / L));
    sqrt(m / L) * (1.0 - sqrt(m / L)) / (1.0 + sqrt(m / L)), 1.0 - sqrt(m / L) + 2.0 * m / L / (1.0 + sqrt(m / L))];

B = [-1.0 / L, -1.0 / L;
    -1.0 / (sqrt(m / L) * L), -1.0 / (sqrt(m / L) * L)];
C = [    (1.0 - sqrt(m / L)) / (1.0 + sqrt(m / L)),2.0 * sqrt(m / L) / (1.0 + sqrt(m / L));   
    sqrt(m / L) * (1.0 - sqrt(m / L)) / (1.0 + sqrt(m / L)),    1.0 - sqrt(m / L) + 2.0 * m / L / (1.0 + sqrt(m / L))];
    
D = [0.0, 0.0;
    -1.0 / (sqrt(m / L) * L), -1.0 / (sqrt(m / L) * L)];

K = ss(A, B, C, D, 1);

%form the system
sys = opt_system(ops, [], K);

% order = [2, 2];
order = {3, 3};

man = opt_analysis(sys);

sol = man.bisect(order);

best_rho = 1 - sqrt(m/L);
[best_rho, sol.rho]


% order = {1, 1};


%% use these IQCs as a warm start in synthesis
iqc = sol.cert.iqc_op;
config = opt_config();
%gradient for f1, proximal for f2
config.syn.D_mask = [0, 0; 1, 1];
man_syn = opt_synthesis(sys, config);
% sol_syn = man_syn.bisect(iqc)
[sh, vh] = man_syn.alternate(3, order, iqc)