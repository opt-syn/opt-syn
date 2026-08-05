%describe the operators 
m = 1; L = 50;
op1 = op_sml(m, L); %gradient of f
op2 = op_pcc();     %indicator function of L1 ball


%run the synthesis routine, use bisection to minimize the convergence rate
sys = opt_system({op1, op2});

man = opt_synthesis(sys); 
sol = man.bisect();

rho = sol.rho; % 0.9050


%% with time-delay dynamics
delay = [1, 0];
network = bridge_channel_delay(delay, delay);
sys = opt_system({op1, op2}, network);
man = opt_synthesis(sys);
sol_delay = man.bisect();

%% analyze performance in the time delay setting
sys_nom = sol.sys;
sys_nom.P = network;
sys_delay = sol_delay.sys;

man_ana_nom = opt_analysis(sys_nom);
man_ana_delay = opt_analysis(sys_delay);
order = {2, 2};
sol_ana_delay = man_ana_delay.bisect(order);

sol_ana_nom = man_ana_nom.bisect(order);

