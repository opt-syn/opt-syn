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

sys = opt_system_switched(ops, network, [], Gring);

reg = regulator_switched(sys);

Pi = cat(2, reg.Pi{:})
config =opt_config();
config.syn.D_mask = 0;
man= opt_synthesis(sys, config);
sol= man.bisect();



%% begin simulation (attempted)
d = 30;
Q = rand_quad(d, m, L);
bstar = randi(101, [d, 1]) - 50;
% bstar = zeros(d, 1);
op1 = op_sim_quad(Q, bstar);

ops_sim = {op1};

sys = sol.sys.export_sim(ops_sim);

T = 50;

sim  = alg_sim(sys, d);
nx = sys.get_alg(1).nx;
sim.sampler.x0 = 4*randn(nx, d);
sim_out = sim.sim(T);

plt = alg_plotter(sim_out);
% plt.plot({'f', 'w', 'res_w', ...
%     'x', 'z', 'mode', }, 10)


% plt.plot({'x', 'w', 'mode', 'z' }, 1)
% plt.plot({'x', 'w', 'mode', 'z' }, 1)
plt.plot({'x', 'mode'}, 1)


%% regulator equations
sim_out_long = sim.sim(20*T);
dstar = -sim_out_long.z(:, :, end);
plt = plt.add_opt_sig(sol.regcl, dstar);

%% new
% plt.plot_4_err(2);    %plot the tracking error
% plt.plot_4_sq_err(3); %plot the  squared norm of the tracking error
plt = plt.add_opt_sig(sol.regcl, dstar);
% plt.plot({'xerr', 'z', 'mode'})
plt.plot_4_err()