%tracking of optimal solution

m = [1, 1, 1, 1];
L = [2, 4, 6, 8];
s = length(m);


omega = pi/8;
Sbeta = blkdiag(1, givens(cos(omega), sin(omega)));
Rbeta = [1, 1, 0];


d = 4; %dimensionality of the problem
%low for illustration


%tracking of the solution
SA = kron(Sbeta, eye(d));
SAy = kron(Rbeta, eye(d));
eta0 = [randi(201, [d, 1]) - 100; 
    (randi(31, [2*d, 1]) + 60).*sign(2*rand(2*d, 1) - 1)];

shift = @(k) SAy * (SA^k) * eta0;

%% assign the quadratics
M = cell(s, 1);
bstar_center = cell(s, 1);
bstar = cell(s, 1);
op_list = cell(s, 1);


for i = 1:s
    M{i} = rand_quad(d, m(i), L(i));
    bstar_center{i} = randi(21, [d, 1]) - 10;    
    bstar{i} = @(k) bstar_center{i} - shift(k);

    op_list{i}= op_sim_quad(M{i}, bstar{i});   
end

%static parameters
b0 = [-0.05; -0.1; -0.05];
b1 = -0.05;
b2 = -0.1;

AK = Sbeta;
BK = b0 * ones(1, s);
CK = ones(s, 1) * Rbeta;

DK = zeros(s);
DK(end, :)= DK(end, :) + b1;
DK(:, 1)= DK(:, 1) + b2;


K = ss(AK, BK, CK, DK, 1);


sys = opt_system(op_list, [], K);

%% execute the algorithm

T = 80;
sim = alg_sim(sys, d);
ssim= sim.sim(T);

%% plot the outputs
plt = alg_plotter(ssim);
plt.plot_6f();
