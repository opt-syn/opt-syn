Nblock = 6;

%symmetry generator/permutation matrix


%define the quadratic


NL = 300;
% NL = 2;
Nblock = 6;
L = logspace(-2, 3, NL);
rho = zeros(NL, Nblock);

parfor i = 1:NL
    for c = 1:Nblock     

        m = 1;
        op1 = op_sml(m, L(i));        

        gamma = 2/(L(i) +m) * (1/c);
        K = ss(eye(c), blkdiag(-gamma , zeros(c-1)), ...
            blkdiag(1, zeros(c-1)), zeros(c), 1);

        M = circshift(eye(c), 1); 

        network = coordinate_descent_primitive(c);
        %form the system

        sys = opt_system_periodic_orbit(op1, network, K, M);

        man = opt_analysis(sys);
        order = [1, 1];
        
        sol = man.bisect(order);
        rho(i, c) = sol.rho;

    end
end

save('ana_coord_sweep_data.mat')

%% simulate and plot
figure(1)
clf
plot(omegalist, rholist_0_0, 'linewidth', 2, 'color', cc(1, :))
plot(omegalist, rholist_1_0, 'linewidth', 2, 'color', cc(2, :))
plot(omegalist, rholist_2_0, 'linewidth', 2, 'color', cc(3, :))
plot([-pi, pi], [1, 1], ':', 'linewidth', 2, 'color', 0.5*[1,1,1])
xlabel('$\omega$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\rho$', 'interpreter', 'latex', 'fontsize', 16)
xlim([min(omegalist), max(omegalist)]);
legend({'Order [0,0]', 'Order [1, 0]', 'Order [2, 0]', '$\rho$=1'}, 'location','northeast', 'interpreter', ...
    'latex', 'fontsize', 16,'location', 'north')
