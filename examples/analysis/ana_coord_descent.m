c=4; %blocksize

%symmetry generator/permutation matrix
M = circshift(eye(c), 1); 

%operator definition
m = 1; L = 5;
% op1 = op_quad(m, L);
op1 = op_sml(m, L);

%gradient descent rule
%coordinate updates
gamma = 2/(L +m) * (1/c); 
K = ss(eye(c), blkdiag(-gamma , zeros(c-1)), ...
    eye(c), zeros(c), 1);

%form the system
sys = opt_system_periodic_orbit(op1, [], K, M);

% %simulate and plot
man = opt_analysis(sys);

order = [2, 1]; %quad: 0.8942 %sml:  0.8942

order = [1, 0];
sol = man.bisect(order) %0.9877
