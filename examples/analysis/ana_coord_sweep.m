Nblock = 6;

%symmetry generator/permutation matrix


%define the quadratic
m = 1; L = 5;
op1 = op_quad(m, L);

NL = 300;
% NL = 2;
Nblock = 6;
L = logspace(-2, 3, NL);
rho = zeros(NL, Nblock);

parfor i = 1:NL
    for c = 1:Nblock     

        gamma = 2/(L(i) +m) * (1/c);
        K = ss(eye(c), blkdiag(-gamma , zeros(c-1)), ...
            blkdiag(1, zeros(c-1)), zeros(c), 1);

        M = circshift(eye(c), 1); 

        network = coordinate_descent_primitive(c);
        %form the system

        sys = opt_system_periodic_orbit(op1, network, K, M);

        man = opt_analysis(sys);
        order = [1, 0];
        
        sol = man.bisect(order);
        rho(i, c) = sol.rho;

    end
end

%% simulate and plot
