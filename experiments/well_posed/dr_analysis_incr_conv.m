%analysis of douglas-rachford algorithm
m= 1;
L = 5;

%different operator classes for op1
op1_sml = op_sml(m, L);
op1_gen = op_gen(1);
op1_gen.monotone = m;
op1_gen.cocoercive = 1/L;

%indicator function for op1
op2 = op_sml(0, inf);

%douglas-rachford
gamma = 0.4;    %stepsizes
lambda = 0.3;
sK = ss([1], [-gamma*lambda, -gamma*lambda], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);


sys_sml = opt_system({op1_gen,  op2}, [], sK);

man_sml = opt_analysis(sys_sml);

sol_sml = man_sml.bisect({0, [2, 0]});
sol_sml_incr = man_sml.bisect({0, [0, 0]});

[sol_sml.rho, sol_sml_incr.rho]