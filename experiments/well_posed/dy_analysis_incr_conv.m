%analysis of douglas-rachford algorithm
m1 = 1;
b1 = 0.1;


op1_gen = op_gen(1);
op1_gen.monotone = m1;
op1_gen.cocoercive = b1;

m2= 2;
L2 = 10;

op2_sml = op_sml(m2, L2);



op3 = op_sml(0, inf);

%douglas-rachford
gamma = 0.3;    %stepsizes
lambda = 0.2;
sK = ss([1], [-gamma*lambda, -gamma*lambda, -gamma*lambda], [1; 1; 1],...
    [-gamma, 0, 0; -gamma, 0, 0; -2*gamma, -gamma, -gamma],1);

sys_sml = opt_system({op1_gen, op2_sml,  op3}, [], sK);

man_sml = opt_analysis(sys_sml);

% sol_iqc_noncaus_bit = man_sml.bisect({3, [3, 3], [3, 3]}); % 0.9729
sol_iqc_noncaus = man_sml.bisect({0, [2, 2], [2, 2]}); %0.9731
sol_incr = man_sml.bisect(0); %0.9731

% sol_iqc_noncaus = man_sml.bisect({[3, 3], 3, [3, 3]}); % 0.9645
% sol_iqc_noncaus = man_sml.bisect({[3, 3], 0, [3, 3]}); %0.9887
% sol_incr = man_sml.bisect(0); %0.9891

[sol_incr.rho, sol_iqc_noncaus.rho, sol_iqc_noncaus_bit.rho]