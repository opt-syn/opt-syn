%sweep of L3
NL = 100;
L_all = logspace(-2, 3, NL);

%define base operators
m = [0, 1, 0];
L = [inf, 1.5, inf];
ops= cell(1, 3);
for i = 1:3
    ops{i} = op_sml(m(i), L(i));    
end

%outer loops
orderlist = {0, 1, [0, 1], [1, 1], 2, [0, 2]};
delaylist = 0:4;
Ndelay = length(delaylist);
Norder = length(orderlist);

rho_best = zeros(NL, 1);
rho_bound = zeros(NL, Norder, Ndelay);
parfor j = 1:NL
    for k = 1:Ndelay
        delay = delaylist(k);
        ops_curr = ops;        
        ops_curr{3}.L = L_all(j);
        
        %davis-yin algorithm
        lambda = 1;        
        gamma = 1/L(2); % or gamma = 1/(2*L(2));
        K = ss([1], [-gamma*lambda, -gamma*lambda, -gamma*lambda], [1; 1; 1], ...
            [-gamma, 0, 0; ...
            -gamma, 0, 0; ...
            -2*gamma, -gamma, -gamma],1);

        %form the system
        network = bridge_channel_delay([0, delay, 0], [0, delay, 0]);
        sys = opt_system(ops_curr, network, K);
    
        %perform analysis at each order
        man = opt_analysis(sys);
        for i = 1:Norder        
            order_curr = repmat({orderlist{i}}, 1, 3);
            sol = man.bisect(order_curr);               
            rho_bound(j, i, k) = sol.rho;    
        end
    end
end


%% plot the result


figure(k)
clf
plot(L_all, squeeze(rho_bound(:, 4, :)), 'linewidth', 2);
hold on
    set(gca, 'Xscale', 'Log')
    ylabel('$\rho$', 'interpreter', 'latex','Fontsize', 14)
    xlabel('$L_3$', 'interpreter', 'latex','Fontsize', 14)
    xlim([L_all(1), L_all(end)])

plot(L_all, ones(size(L_all)), ':', 'LineWidth', 2, 'color', 0.5*[1,1,1])