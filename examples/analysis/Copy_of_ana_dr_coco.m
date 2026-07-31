%cocoercive plus strongly monotone
mu = 1;
beta = 1.5;


SUBDIFF = false;

if SUBDIFF
    op1 = op_sml(mu, inf);
    op2 = op_sml(0, 1/beta);
else    
    op1 = op_gen();
    op1.monotone = mu;
    
    op2 = op_gen();
    op2.cocoercive= beta;
end

if beta^2 + mu^beta + beta - mu <= 0
    rhobest = 1 - beta/(1+beta);
elseif mu^beta - mu - beta >= 1
    rhobest = 1 - (1+mu*beta)/((mu+1)*(beta+1));
elseif mu^2 + beta*mu + mu - beta <= 0
    rhobest = 1 - mu/(mu+1);
else
    rhobest = 1/2 * (beta+mu)/sqrt(mu*beta*(beta + mu + 1));
end
    



%douglas rachford
gamma = 1;
lambda = 1;

K = ss([1], [-lambda*gamma, -lambda*gamma], [1; 1], [-gamma, 0; -2*gamma, -gamma],1);

sys = opt_system({op1, op2}, [], K);

man = opt_analysis(sys);
order = {0, 1};

% specs = spec_stability(0.3);
% sol = man.solve_single(order, specs)
sol_ana = man.bisect(order)

% man.bisect(order)