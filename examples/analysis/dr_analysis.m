
%define the quadratic
m = 1; L = 10;
op1 = op_sml(m, L);
op2 = op_pcc();

%projected pseudogradient descent
gamma = 0.4;
K = ss([1], [-gamma, -gamma], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);

sys = opt_system({op1, op2}, [], K);

man = opt_analysis(sys);
order = {[1, 1], [1, 1]};
man.bisect(order)