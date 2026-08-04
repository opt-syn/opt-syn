%analysis of douglas-rachford algorithm


m= 1;
L = 10;

%different operator classes for op1
op1_gen = op_gen({'monotone', 1, 'cocoercive', 1/L});
op1_sml = op_sml(m, L);
op1_quad = op_quad(m, L);

%indicator function for op1
op2 = op_sml(0, inf, 1);


%douglas-rachford
gamma = 0.4;    %stepsizes
lambda = 0.25;
sK = ss([1], [-gamma*lambda, -gamma*lambda], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);


sys_gen = opt_system({op1_gen, op2}, [], sK);
sys_sml = opt_system({op1_sml,  op2}, [], sK);
sys_quad = opt_system({op1_quad, op2}, [], sK);

man_gen = opt_analysis(sys_gen);