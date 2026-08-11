%define the operators
m = 1;L = 10;
op1= op_pcc(); %pcc goes first, to not take information from the noisy gradient
op2= op_sml(m, L);
ops = {op1, op2};

%define the controller (unconstrained and constrained/composite)
gamma = 2/(L + m);
K = ss(1, [-gamma, -gamma], [1; 1], [-gamma, 0; -gamma, 0],1);

%gain from (error in w2) to (z2 - z*)
network = bridge_pass_through(2);
[network, iwp] = network.add_oracle_input(2, []);
[network, izp] = network.perf_output_z(2);



%define specification
GAIN = 10;
spec = spec_h2(GAIN, 1, iwp, izp);
spec.target = true;
specs = {spec};

%form the systems
sys_clean = opt_system(ops, [], K);
sys_noisy = opt_system(ops, network, K);

%pose the analysis problem
order = [1,1];

man_noisy= opt_analysis(sys_noisy);
sol_noisy = man_noisy.solve_single(order, specs);

man_clean = opt_analysis(sys_clean);
sol_clean = man_clean.bisect(order);


l2gain = sqrt(sol_noisy.objective);
rho = sol_clean.rho;