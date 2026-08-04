%time-varying delays
rng(40, 'twister');
d = 64;
delay_max = 2;


%% define the operators
%generate the test quadratic
m = 1; L = 1.5;
ops = {op_sml(m, L)};


%% form the time-varying delay network

%delay for oracle 1
% [Pprim, Gcon] = delay_primitives(0, 0:delay_max, 0,  -1:1);
z = tf('z');
alpha = 0.9;
network = cell(delay_max+1, 1);
for i = 1:delay_max+1
    network{i} = genplant([0, 1; (z/(z-alpha^(i-1)))^2, 0]);
    network{i}.nw = 1; network{i}.ny = 1;
    network{i}.nz = 1; network{i}.nu = 1;
end
Gcon = toeplitz([1; 1; zeros(delay_max-1, 1)], [1; 1; zeros(delay_max-1, 1)]);
Gsnap = delay_snap_graph(delay_max, 0);

Nss = length(Gsnap);

network = genplant_poly(network);

%gradient descent

%form the systems
sys_per = opt_system_periodic(ops, network, []);
sys_snap = opt_system_switched(ops, network, [], Gsnap);
sys_cont = opt_system_switched(ops, network, [], Gcon);

%% synthesize
config =opt_config();
config.syn.D_mask = 0;
man_per = opt_synthesis(sys_per, config);
man_snap = opt_synthesis(sys_snap, config);
man_cont = opt_synthesis(sys_cont, config);

sol_per = man_per.bisect();
sol_snap = man_snap.bisect();
sol_cont = man_cont.bisect();


%% begin simulation (attempted)
Q = rand_quad(d, m, L);
bstar = randi(101, [d, 1]) - 50;
op1 = op_sim_quad(Q, bstar);

ops_sim = {op1};

sys_per = sol_per.sys.export_sim(ops_sim);
sys_snap = sol_snap.sys.export_sim(ops_sim);
sys_cont = sol_cont.sys.export_sim(ops_sim);

T = 150;
sim_per  = alg_sim(sys_per, d);
sim_snap = alg_sim(sys_snap, d);
sim_cont = alg_sim(sys_cont, d);


sim_out_per  = sim_per.sim(T);
sim_out_snap = sim_snap.sim(T);
sim_out_cont = sim_cont.sim(T);

%% plot the signals
plt_per = alg_plotter(sim_out_per);
plt_per.plot({'f', 'w', 'res_w', ...
    'x', 'z', 'delay', }, 10)

plt_snap = alg_plotter(sim_out_snap);
plt_snap.plot({'f', 'w', 'res_w', ...
    'x', 'z', 'delay', }, 12)


plt_cont = alg_plotter(sim_out_cont);
plt_cont.plot({'f', 'w', 'res_w', ...
    'x', 'z', 'delay', }, 11)
