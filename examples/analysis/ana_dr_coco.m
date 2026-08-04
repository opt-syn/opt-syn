%cocoercive plus strongly monotone
mu = 1; beta = 1.5;

%sweep beta
Nbeta = 100;
beta_sweep = logspace(-2, 2, Nbeta);

%description of operators
op1 = op_gen();
op1.monotone = mu;

op2 = op_gen();

%douglas rachford
gamma = 1; lambda = 1;

K = ss([1],    [-lambda*gamma, -lambda*gamma], ...
       [1; 1], [-gamma, 0;  -2*gamma, -gamma], 1);

%form the system
sys = opt_system({op1, op2}, [], K);

%solve the problem
order = {1, 1};
man = opt_analysis(sys);

rho_best_dr = zeros(Nbeta, 1);
rho = zeros(Nbeta, 1);
%set and solve
for i = 1:Nbeta
    %assign the sweep
    man.sys.op{2}.cocoercive = beta_sweep(i);    

    %compute bound
    sol = man.bisect(order);
    rho(i) = sol.rho;   

    %find the true rate
    rho_best_dr(i) = best_dr_rate(mu, beta_sweep(i));     
end



%% plot the output
cc = linspecer(4);
figure(3)
clf
tiledlayout(2, 1)
nexttile
hold on
plot(beta_sweep, rho,  'LineWidth', 2, 'color', cc(1, :))
plot(beta_sweep, rho_best_dr, ':', 'LineWidth', 2, 'color', cc(2, :))
set(gca, 'Xscale', 'log')
title('Contraction Rates', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\rho$', 'interpreter', 'latex')
xlabel('$\beta$', 'interpreter', 'latex')
legend({'upper-bound', 'best rate'}, 'location','northeast', 'interpreter', 'latex')

nexttile
plot(beta_sweep, rho - rho_best_dr, 'LineWidth', 2, 'color', cc(3, :))
set(gca, 'Xscale', 'log')
xlabel('$\beta$', 'interpreter', 'latex')
ylabel('$\rho_{bound} - \rho_{best}$', 'interpreter', 'latex')
title('Error in Estimates', 'interpreter', 'latex', 'FontSize', 20)





%% compute best DR rate
function rhobest = best_dr_rate(mu, beta)
    %optimal convergence rate for strong monotone + cocoercive
    if beta^2 + mu*beta + beta - mu <= 0
        rhobest = abs(1 - beta/(1+beta));
    elseif mu^beta - mu - beta >= 1
        rhobest = abs(1 - (1+mu*beta)/((mu+1)*(beta+1)));
    elseif mu^2 + beta*mu + mu - beta <= 0
        rhobest = abs(1 - mu/(mu+1));
    else
        rhobest = 1/2 * (beta+mu)/sqrt(mu*beta*(beta + mu + 1));
    end
end
