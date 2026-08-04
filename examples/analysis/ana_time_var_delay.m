%time-varying delays
rng(40, 'twister');
d = 64;
delay_max = 3;

%% define the operators
m = 1; L = 1.5;
op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};

%% form the time-varying delay network
%delay for oracle 1
[Pprim, Gcon] = delay_primitives(0, 0:delay_max, -1:1, 0);
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

sys_snap = opt_system_switched(ops, network, K, Gsnap);
sys_cont = opt_system_switched(ops, network, K, Gcon);
sys_per = opt_system_periodic(ops, network, K);

sys_lti = sys_per.periodic_lift();

%% perform analysis
order = {1, 1};

man_lti = opt_analysis(sys_lti);
sol_lti = man_lti.bisect(order);

man_per = opt_analysis(sys_per);
sol_per = man_per.bisect(order);

man_snap = opt_analysis(sys_snap);
sol_snap = man_snap.bisect(order);

man_cont = opt_analysis(sys_cont);
sol_cont = man_cont.bisect(order);

[sol_lti.rho, sol_per.rho, sol_snap.rho, sol_cont.rho]
%0.8000    0.9460    0.9500    0.9500

