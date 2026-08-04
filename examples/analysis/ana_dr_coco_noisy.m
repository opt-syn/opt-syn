%cocoercive plus strongly monotone
mu = 1; beta = 1.5;

%sweep beta
% Nbeta = 100;
Nbeta = 3;
beta_sweep = logspace(-2, 2, Nbeta);

%description of operators
op1 = op_gen();
op1.monotone = mu;

op2 = op_gen();

%douglas rachford
gamma = 1; lambda = 1;

K = ss([1],    [-lambda*gamma, -lambda*gamma], ...
    [1; 1], [-gamma, 0;  -2*gamma, -gamma], 1);

%create the noise and readouts
network = bridge_pass_through(2);
[network, iwp] = network.add_oracle_input([], [1, 2]);
[network, izp] = network.perf_output_opt();

%form the system
sys = opt_system({op1, op2}, network, K);


%solve the problem
order = {1, 1};
man = opt_analysis(sys);

rho_best_dr = zeros(Nbeta, 1);
rho = zeros(Nbeta, 1);
%set and solve



spec = spec_e2e(100, iwp, izp);
spec.target = true;
GAIN = zeros(Nbeta, 1);
for i = 1:Nbeta
    %assign the sweep
    man.sys.op{2}.cocoercive = beta_sweep(i);    

    %compute bound
    sol = man.solve_single(order, {spec});
    GAIN(i) = sol.objective;   

    %find the true rate
end
