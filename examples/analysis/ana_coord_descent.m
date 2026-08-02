c=3; %blocksize

% c=6; %blocksize

%symmetry generator/permutation matrix
M = circshift(eye(c), 1); 

%gradient descent rule
gamma = 0.1;
K = ss(eye(c), -gamma * eye(c), eye(c), zeros(c), 1);
% K = ss(eye(c), blkdiag(-gamma , zeros(c-1)), ...
%     blkdiag(1, zeros(c-1)), zeros(c), 1);


%define the quadratic
m = 1; L = 5;
% op1 = op_sml(m, L);
op1 = op_quad(m, L);

%define the network
%implements coordinate descent
[Pprim] = coordinate_descent_primitives(c);
network = Pprim{1};

%form the system
sys = opt_system_periodic_orbit(op1, network, K, M);

%regulator equation
% reg = regulator_periodic_orbit(sys);
% rcl = reg.check_regulator()

% %simulate and plot
man = opt_analysis(sys);
order = [1, 0];
%
sol = man.bisect(order) %0.9828