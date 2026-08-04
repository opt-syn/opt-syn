%theorem 2 requires that lambda is in (0, 2-alpha(m1+L1)/2)
%and alpha > 0. means that alpha in (0, 2).

% NL = 10;

% L_all = logspace(-2, 3, NL);
L_all = 3;
NL = 1;

m = [0, 1, 0];
L = [inf, 1.5, inf];
ops= cell(1, 3);
for i = 1:3
    ops{i} = op_sml(m(i), L(i));    
end
% ops{1}.ERGODIC = true;

orderlist = {0, 1, [1, 1]};
Norder = length(orderlist);

rho_best = zeros(NL, 1);
rho_bound = zeros(NL, Norder);
for j = 1:NL
    L(3) = L_all(j);
    ops{3}.L = L(3);
    [rho_best(j), K] = best_dy(m, L);

    
    sys = opt_system(ops, [], K);

    man = opt_analysis(sys);
    for k = 1:Norder
        order_curr = repmat({orderlist{k}}, 1, 3);
        sol = man.bisect(order_curr);               
        rho_bound(j, k) = sol.rho;    
    end
end



%% ergodic performance
[perf_erg, sys_erg] = spec_ergodic(sys);
specs = {perf_erg};
sys_erg.op{3}.L = inf;
sys_erg.op{1}.ERGODIC = true;

man_erg = opt_analysis(sys_erg);

% sol_erg = man.solve_single(order, specs);



%% plot the result
figure(5)
clf
hold on
plot(L_all, rho_bound, 'linewidth', 2)
plot(L_all, rho_best, 'linewidth', 2)
ylim([0, 1])
set(gca, 'Xscale', 'Log')
xlim([L_all(2), L_all(end)])
ylim([min(rho_bound), 1])


%% best bound
function [rho, K] = best_dy(m, L)
    %compute the optimal davis-yin rate and controller
    
    
    
    Nalpha = 10000;
    alphalist = linspace(0, 2, Nalpha+2);
    alphalist = alphalist(2:end-1);
    lambdalist = zeros(Nalpha, 1);
    
    rholist = zeros(Nalpha, 1);
    for i = Nalpha:-1:1
        alpha = alphalist(i);
    
    
        nu_g_1 = (2*m(3) + m(2))/(1+alpha * m(3))^2;
        nu_g_2 = (2*L(3) + m(2))/(1+alpha * L(3))^2;
    
        nu_g = min(nu_g_1, nu_g_2);
        theta = 2/(4-alpha * (m(2) + L(2)));
    
        gap = sqrt(theta * (theta - alpha*nu_g)) - theta;
    
    
        if imag(gap)==0        
            lambda = (1-1e-5)*(2 - alpha*(m(2) + L(2))/2) * (gap< 0);
            rholist(i) = sqrt(1 + lambda*(gap));
        else
            rholist(i) = NaN;
        end
    
        lambdalist(i) = lambda;
    
    
    
    end
    
    
  
    [rho, midx] = min(rholist, [], 'all');
    lambda = lambdalist(midx);
    alpha = alphalist(midx);

    

    K = ss([1], [-alpha*lambda, -alpha*lambda, -alpha*lambda], [1; 1; 1], ...
        [-alpha, 0, 0; -alpha, 0, 0;-2*alpha, -alpha, -alpha],1);

end