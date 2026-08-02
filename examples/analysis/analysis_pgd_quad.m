%analysis of douglas-rachford algorithm
m= 1;
L = 5;

%different operator classes for op1
op1_sml = op_sml(m, L);
op1_quad = op_quad(m, L);

%indicator function for op1
op2 = op_sml(0, inf, 1);


%projected gradient descent
gamma = 2/(L+m);    
sK = ss([1], [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);

network = [];
network = bridge_channel_delay([1, 0], [1, 0]);

sys_sml = opt_system({op1_sml,  op2}, network,  sK);
sys_quad = opt_system({op1_quad, op2}, network, sK);

order = [2, 2];

man_sml = opt_analysis(sys_sml);
man_quad = opt_analysis(sys_quad);

sol_sml = man_sml.bisect(order);
sol_quad = man_quad.bisect(order);


%% analyze solution

P_sml = op1_sml.dhd_lift(order, sol_sml.vars.op{1}, sol_sml.cert.iqc_op{1});
P_quad = op1_quad.dhd_lift(order, sol_quad.vars.op{1}, sol_quad.cert.iqc_op{1});


