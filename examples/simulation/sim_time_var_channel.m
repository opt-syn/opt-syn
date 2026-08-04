%time-varying delays
rng(40, 'twister');
d = 64;
delay_max = 3;


%% define the operators
%generate the test quadratic
m = 1; L = 1.5;
Q = rand_quad(d, m, L);
bstar = randi(101, [d, 1]) - 50;
op1 = op_sim_quad(Q, bstar);

ops = {op1};

%% form the time-varying delay network

%delay for oracle 1
% [Pprim, Gcon] = delay_primitives(0, 0:delay_max, 0,  -1:1);
z = tf('z');
alpha = 0.9;
network = cell(delay_max+1, 1);
for i = 1:delay_max+1
    network{i} = genplant([0, z/(z-alpha^(i-1)); 1, 0]);
    network{i}.nw = 1; network{i}.ny = 1;
    network{i}.nz = 1; network{i}.nu = 1;
end
Gcon = toeplitz([1; 1; zeros(delay_max-1, 1)], [1; 1; zeros(delay_max-1, 1)]);
Gsnap = delay_snap_graph(delay_max, 0);

Nss = length(Gsnap);

network = genplant_poly(network);

%gradient descent
gamma = 0.05;
K = ss([1], [-gamma], [1],0,1);

%form the systems
sys_per = opt_system_periodic(ops, network, K);
sys_snap = opt_system_switched(ops, network, K, Gsnap);
sys_cont = opt_system_switched(ops, network, K, Gcon);


%% begin simulation (attempted)
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
