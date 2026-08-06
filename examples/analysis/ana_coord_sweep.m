Nblock = 6;



%define the quadratic
NL = 300;
Nblock = 6;
L = logspace(-2, 3, NL);
rho = zeros(NL, Nblock);

parfor i = 1:NL
    for c = 1:Nblock     

        m = 1;
        op1 = op_sml(m, L(i));        

        %gradient descent scheme
        gamma = 2/(L(i) +m) * (1/c);
        K = ss(eye(c), blkdiag(-gamma , zeros(c-1)), ...
            blkdiag(1, zeros(c-1)), zeros(c), 1);

        %symmetry generator/permutation matrix

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
hold on
plot(L(121:end), rho(121:end, :), 'linewidth', 2)
% % plot(omegalist, rholist_1_0, 'linewidth', 2, 'color', cc(2, :))
% plot(omegalist, rholist_2_0, 'linewidth', 2, 'color', cc(3, :))
plot([L(1), L(end)], [1, 1], ':', 'linewidth', 2, 'color', 0.5*[1,1,1])
xlabel('$c$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\rho$', 'interpreter', 'latex', 'fontsize', 16)
xlim([1, 100]);
ylim([0, 1.2])
Lname = cell(Nblock+1, 1);
for i = 1:Nblock
    Lname{i} = ['$c=', num2str(i), '$'];
end
Lname{end} = '$\rho=1$';
legend(Lname, 'location','northeast', 'interpreter', ...
    'latex', 'fontsize', 16,'location', 'southeast')

set(gca, 'xscale', 'log')