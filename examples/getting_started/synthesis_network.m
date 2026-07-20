% The set of considered oracles
m = 1;
L = 10;
op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};

DELAY = [2, 0];
network_delay = bridge_channel_delay(DELAY, DELAY);


%form the system and the analysis manager
sys_delay = opt_system(ops, network_delay);
man_delay = opt_synthesis(sys_delay);
sol_best_delay = man_delay.bisect();


%% plot an execution
d = 100;
M = rand_quad(d, m, L);
zstar = 100*randn(d);
op1_sim = op_sim_quad(M, zstar);
op2_sim = op_sim_box(30);

sys_sim = sys_delay;
sys_sim.op = {op1_sim, op2_sim};