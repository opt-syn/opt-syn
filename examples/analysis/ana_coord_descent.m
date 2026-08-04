c=3; %blocksize

%symmetry generator/permutation matrix
M = circshift(eye(c), 1); 

%operator definition
m = 1; L = 5;
op1 = op_quad(m, L);


%gradient descent rule
gamma = 2/(L +m) * (1/c); %0.9685
K = ss(eye(c), blkdiag(-gamma , zeros(c-1)), ...
    blkdiag(1, zeros(c-1)), zeros(c), 1);

%define the network
network = coordinate_descent_primitive(c);


%form the system
sys = opt_system_periodic_orbit(op1, network, K, M);

% %simulate and plot
man = opt_analysis(sys);

order = [1, 0];
sol = man.bisect(order) %0.9712