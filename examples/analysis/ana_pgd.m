%analysis of douglas-rachford algorithm
m= 1;
L = 10;

%different operator classes for op1
op1_sml = op_sml(m, L);

%indicator function for op1
op2 = op_sml(0, inf);


%projected gradient descent
gamma = 2/(L+m); 
% gamma = 1/L;
sK = ss([1], [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);


sys_sml = opt_system({op1_sml,  op2}, [],  sK);

order = {[0, 0], [0, 0]};

config = opt_config();
config.recovery.blocks = false;
man_sml = opt_analysis(sys_sml, config);


sol_sml = man_sml.bisect(order);
sol = sol_sml;
rho_sml = sol_sml.rho
%% analyze solution

P_sml = op1_sml.dhd_lift(order{1}, sol_sml.vars.op{1}, sol_sml.cert.iqc_op{1});
H_sml = P_sml + P_sml' 