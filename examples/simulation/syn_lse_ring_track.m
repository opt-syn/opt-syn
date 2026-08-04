%generate the subsystems
z = tf('z', 1);

%original
D11_base = {0, 0, 0, 0};
D12_base = {1/(z-0.2), 1/(z-0.9), 1/(z+0.5), 1/(z+0.2)};
D21_base = {-1/z, -0.5*z/(z-0.2), 0.5*z/(z+0.3), 1/(z-1.2)};
D22_base = {0, 3, 1, 2};

n = struct('nw', 1, 'nz', 1, 'ny', 1, 'nu', 1);

Plist = cell(4, 1);
for i = 1:4
    Plist{i} = genplant(ss([D11_base{i}, D12_base{i}; D21_base{i}, D22_base{i}]), n);
end
network = genplant_poly(Plist);

Gring = [1, 1, 0, 0;
    0, 1, 1, 0;
    0, 0, 1, 1;
    1, 0, 0, 1];

m = 1; L = 1.5;
ops = {op_sml(m, L)};

%tracking
theta = pi/8;    
% theta = pi/30;
Sbeta = blkdiag(1, givens(cos(theta), sin(theta)));
Rbeta = [1, 1, 0];
tracking = struct('Sbeta', Sbeta, 'Rbeta', Rbeta);

sys = opt_system_switched(ops, network, [], Gring);

sys.tracking = tracking;

reg = regulator_switched(sys);

Pi = cat(2, reg.Pi{:})
config =opt_config();
config.syn.D_mask = 0;
man= opt_synthesis(sys, config);
sol= man.bisect();



%% begin simulation (attempted)
d = 40;
Q = rand_quad(d, m, L);

%shift
SA = kron(Sbeta, eye(d));
SAy = kron(Rbeta, eye(d));
dstar0 = [randi(201, [d, 1]) - 100; (randi(31, [2*d, 1]) + 60).*sign(2*rand(2*d, 1) - 1)];

shift = @(k) SAy * (SA^k) * dstar0;

bstar = @(k, param) randi(101, [d, 1]) - 50 + shift(k);


op1 = op_sim_quad(Q, bstar);

ops_sim = {op1};

sys = sol.sys.export_sim(ops_sim);

T = 50;

sim  = alg_sim(sys, d);
sim_out = sim.sim(T);
plt = alg_plotter(sim_out);
plt.plot({'f', 'w', 'res_w', ...
    'x', 'z', 'mode', }, 10)

