%define the operators
m = 1;L = 10;
op1= op_sml(m, L);

%define the controller (unconstrained and constrained/composite)
K_uncon = tmm(m, L);
 

%add noise to the network

%gain from (error in w) to (z - z*)
network = bridge_pass_through(1);
[network, iwp] = network.add_oracle_input(1, []);
[network, izp] = network.perf_output_z(1);



%define specification
GAIN = 10;
spec = spec_h2(GAIN, 1, iwp, izp);
spec.target = true;
specs = {spec};

%form the system
sys_uncon = opt_system(op1, network, K_uncon);
 

%pose the analysis problem
order = [1,1];

man_uncon = opt_analysis(sys_uncon);
 
%solve the analysis problem
sol_uncon = man_uncon.solve_single(order, specs)
% sol_uncon = man_uncon.bisect(order);
 
% best_rho = 1 - sqrt(m/L);
% [best_rho, sol_uncon.rho, sol_con.rho]
% 

%% functions to generate algorithms
function K = tmm(m, L)
%triple momentum algorithm
rho=1-1/sqrt(L/m);
al=(1+rho)/L;
be=rho^2/(2-rho);
ga=rho^2/((1+rho)*(2-rho));
A=[1+be -be;1 0];B=[-al;0];C=[1+ga -ga];
K = ss(A,B,C,0,1);
end
 