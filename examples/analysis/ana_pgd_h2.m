%define the operators
m = 1;L = 10;
op1= op_pcc(); %pcc goes first, to not take information from the noisy gradient
op2= op_sml(m, L);
ops = {op1, op2};

%define the PGD
gamma = 2/(L + m);
K = ss(1, [-gamma, -gamma], [1; 1], [-gamma, 0; -gamma, 0],1);

%gain from (error in w2) to (z2 - z*)
network = bridge_pass_through(2); %start with no network dynamics
[network, iwp] = network.add_oracle_input(2, []); %add wp
[network, izp] = network.perf_output_z(2); %add zp

%define specification
GAIN = 10;
spec_stoch = spec_h2(GAIN, 1, iwp, izp);
spec_stoch.target = true;


%form the systems
sys_clean = opt_system(ops, [], K); %for convergence rate
sys_noisy = opt_system(ops, network, K); %for stochastic sensitivity

%pose the analysis problem
order = {1, 1};

%only convergence
man_noisy= opt_analysis(sys_noisy);
sol_noisy = man_noisy.solve_single(order, spec_stoch);

%only sensitivity
man_clean = opt_analysis(sys_clean);
sol_clean = man_clean.bisect(order);

%both at once: minimize rho with inner objective h2 gain
specs_joint = {spec_stability(1), spec_stoch};
man_joint = opt_analysis(sys_noisy);
sol_joint = man_joint.bisect(order, specs_joint);


h2gain = sqrt(sol_noisy.objective); %  0.3162
rho = sol_clean.rho; % 0.8182

h2gain_joint = sqrt(sol_joint.objective); % 1.0389
rho_joint = sol_joint.rho; % 0.8182