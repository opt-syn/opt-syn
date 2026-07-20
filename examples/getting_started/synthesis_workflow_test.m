% The set of considered oracles
m = 1;
L = 10;
op1 = op_sml(m, L);
op2 = op_pcc();
ops = {op1, op2};

%form the system and the analysis manager
sys = opt_system(ops, [], K);
man = opt_synthesis(sys);
sol_best = man.bisect();


%% gradient
config_grad = opt_config();
config_grad.syn.D_mask = [0, 0; 1, 1];
man_grad = opt_synthesis(sys, config_grad);
sol_grad = man_grad.bisect();

%% alternation
man_alt = man;
b_opts = bisect_opts;
b_opts.Niter = 3;
sol_alt = man_alt.alternate([], {1, 1}, [], b_opts);
% %% delay
% % We introduce a time delay of 2 steps before and after the first oracle evaluation
% DELAY = [2, 0];
% network_delay = bridge_channel_delay(DELAY, DELAY);
% 
% % order_delay = {[2, 2], [2, 2]};
% order_delay = {1, 1};
% sys_delay = opt_system(ops, network_delay, K);
% man_delay = opt_analysis(sys_delay);
% [sol_best_delay, v_range_delay] = man_delay.bisect(order_delay);