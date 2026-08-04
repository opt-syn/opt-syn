%time-varying delays
rng(40, 'twister');
d = 64;
delay_max = 2;

m = 1;
L = 1.5;


%generate the test quadratic
Q = rand_quad(d, m, L);
bstar = randi(101, [d, 1]) - 50;


BOX = 30;

op1 = op_sim_quad(Q, bstar);
op2 = op_sim_box(BOX);
ops = {op1, op2};


%% form the network

[Pprim, Gcon] = delay_primitives(0:delay_max, 0:delay_max, -1:0:1, -1:0:1);

% Gsnap = delay_snap_graph(delay_max, 0);
Gsnap = delay_snap_graph(delay_max, delay_max);
Nss = length(Gsnap);

P0 = bridge_pass_through(1);
P_chain = cell(Nss, 1);
for i = 1:Nss    
    P_chain{i} = blkdiag(Pprim{i}, P0);
end

network = genplant_poly(P_chain);

gamma = 0.4;
lambda = 0.25;
K = ss([1], [-gamma*lambda, -gamma*lambda], [1; 1], [0, 0; -gamma, -gamma],1);

sys_per = opt_system_periodic(ops, network, K);
sys_snap = opt_system_switched(ops, network, K, Gsnap);

%% begin simulation (attempted)
T = 150;
sim = alg_sim(sys_snap, d);

ssim= sim.sim(T);



%% plot the signal
plt = alg_plotter(ssim);
plt.plot({'x', 'z', 'res_z', 'mode', 'w', 'res_w'}, 13)
