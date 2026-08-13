%time-varying delays
rng(40, 'twister');
d = 20;
delay_max = 3;


%% define the operators
%generate the test quadratic
m = 1; L = 10;
Q = rand_quad(d, m, L);
bstar = randi(101, [d, 1]) - 50;
op1 = op_sim_quad(Q, bstar);

ops = {op1};

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
    P_chain{i} = Pprim{i};
end

network = genplant_poly(P_chain);

%projected gradient descent
gamma = 2/(L+m)*0.5;
K = ss([1], [-gamma], [1], [0],1);


% Pl = cellfun(@(s) lft(loop, lft(s.P, K), 1, 1), P_chain, UniformOutput=false);
% ePl = cellfun(@eig, Pl, UniformOutput=false);
% ePl = cat(2, ePl{:})
% AP = cellfun(@(s) s.A, Pl, 'UniformOutput', false);
% [bounds, P, info] = jsr_conic_ellipsoid(AP)
% max(abs(ePl), [], 1)

%form the systems
sys_none = opt_system(ops, [], K);
sys_3 = opt_system(ops, Pprim{end}, K);
sys_per = opt_system_periodic(ops, network, K);
sys_snap = opt_system_switched(ops, network, K, Gsnap);
sys_cont = opt_system_switched(ops, network, K, Gcon);


%% begin simulation
T = 100;
sim_none = alg_sim(sys_none, d);
sim_out_none = sim_none.sim(T);
plt_per = alg_plotter(sim_out_none);
plt_per.plot({'z', 'res_w'}, 10)


sim_3 = alg_sim(sys_3, d);
sim_out_3 = sim_3.sim(T);
plt_per = alg_plotter(sim_out_3);
plt_per.plot({'z', 'res_w'}, 11)



sim_per  = alg_sim(sys_per, d);
sim_snap = alg_sim(sys_snap, d);
sim_cont = alg_sim(sys_cont, d);


sim_out_per  = sim_per.sim(T);
sim_out_snap = sim_snap.sim(T);
sim_out_cont = sim_cont.sim(T);

%% plot the signals
plt_per = alg_plotter(sim_out_per);
plt_per.plot({'delay', 'z', 'res_w'}, 13)

plt_snap = alg_plotter(sim_out_snap);
plt_snap.plot({'delay', 'z', 'res_w'}, 14)


plt_cont = alg_plotter(sim_out_cont);
plt_cont.plot({'delay', 'z', 'res_w'}, 15)
