
%define the quadratic
m = 1; L = 10;
op1 = op_sml(m, L);
op2 = op_pcc();

%douglas-rachford
% gamma = (sqrt(2)-1)/L;    %stepsizes
% lambda = (1-gamma*L)/(1+gamma*L);
gamma = 0.4;
lambda = 1;
K = ss([1], [-gamma*lambda, -gamma*lambda], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);

sys = opt_system({op1, op2}, [], K);

man = opt_analysis(sys);
order = {[1, 1], [1, 1]};
man.bisect(order)