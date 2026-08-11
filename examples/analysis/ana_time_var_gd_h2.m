%time-varying delays
rng(40, 'twister');
d = 64;
delay_max = 1;

%% define the operators
m = 1; L = 1.5;
op1 = op_sml(1, L);
 
ops = {op1};

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
    Pcurr =  Pprim{i} ;
    [Pcurr, iwp] = Pcurr.add_oracle_input(1, []);
    [Pcurr, izp] = Pcurr.perf_output_z(1);

    P_chain{i} = Pcurr;
end

network = genplant_poly(P_chain);

%projected gradient descent
% gamma = 0.05;
gamma = 2/(L+m);
K = ss([1], -gamma, 1, 0,1);

sys_snap = opt_system_switched(ops, network, K, Gsnap);
sys_cont = opt_system_switched(ops, network, K, Gcon);
sys_per = opt_system_periodic(ops, network, K);

sys_lti = sys_per.periodic_lift();

%% perform analysis
order = {[1, 1], [1, 1]};

spec = spec_h2(10, 1, 1, 1);
spec.target = true;
specs = {spec};


% man_lti = opt_analysis(sys_lti);
% sol_lti = man_lti.solve_single(order, specs);

man_per = opt_analysis(sys_per);
sol_per = man_per.solve_single(order, specs);

% man_snap = opt_analysis(sys_snap);
% sol_snap = man_snap.solve_single(order, specs);
% 
% man_cont = opt_analysis(sys_cont);
% sol_cont = man_cont.solve_single(order, specs);

[sol_per.objective, sol_snap.objective, sol_cont.objective]


