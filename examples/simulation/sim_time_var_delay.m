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

%indicator function 
BOX = 30;
op2 = op_sim_box(BOX);

ops = {op1, op2};

%% form the time-varying delay network

%delay for oracle 1
[Pprim, Gcon] = delay_primitives(0, 0:delay_max, 0,  -1:1);
Gsnap = delay_snap_graph(delay_max, 0);

%no delay for oracle 2
P0 = bridge_pass_through(1);

%build the network by block diagonalization
Nss = length(Gsnap);
P_chain = cell(Nss, 1);
for i = 1:Nss    
    P_chain{i} = blkdiag(Pprim{i}, P0);
end

network = genplant_poly(P_chain);

%projected gradient descent
gamma = 0.05;
K = ss([1], [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);

%form the systems
sys_per = opt_system_periodic(ops, network, K);
sys_snap = opt_system_switched(ops, network, K, Gsnap);
sys_cont = opt_system_switched(ops, network, K, Gcon);


%% begin simulation
T = 150;
sim_per  = alg_sim(sys_per, d);
sim_snap = alg_sim(sys_snap, d);
sim_cont = alg_sim(sys_cont, d);


sim_out_per  = sim_per.sim(T);
sim_out_snap = sim_snap.sim(T);
sim_out_cont = sim_cont.sim(T);

%% plot the signals
plt_per = alg_plotter(sim_out_per);
plt_per.plot({'x', 'z', 'res_z', 'delay', 'w', 'res_w'}, 10)

plt_snap = alg_plotter(sim_out_snap);
plt_snap.plot({'x', 'z', 'res_z', 'delay', 'w', 'res_w'}, 12)


plt_cont = alg_plotter(sim_out_cont);
plt_cont.plot({'x', 'z', 'res_z', 'delay', 'w', 'res_w'}, 11)
