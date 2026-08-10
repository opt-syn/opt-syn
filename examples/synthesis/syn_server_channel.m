%form the operators
m = [0, 1, -2, 1, 1, 0];
L = [5, 2, 1, 6, 1, inf];
s = length(m);

ops= cell(s, 1);
for i = 1:s
    ops{i} = op_sml(m(i), L(i));
end

%describe the network
channel = bridge_pass_through(s);
z = tf('z', 1);
channel.P(5, s+2) = z^(-1);
channel.P(s+2, 5) = z^(-1);

channel.P(s+4, 4) = z^(-2);
channel.P(s+3, 5) = z^(-1);

channel.P(1, s+1) =  0.5/(z-0.5);
channel.P(s+1, 1) = 0.5/(z-0.5);

channel.P(3, s+3) =  -1/(z+0.4);
channel.P(s+3, 3) = -1/(z+0.4);

%form the systems
sys = opt_system(ops);
sys_channel = opt_system(ops, channel);


%information structure for the channel case
config = opt_config();
config.syn.D_mask = [1, 0, 1, 1, 0, 0;
    1, 1, 1, 1, 0, 0;
    1, 1, 1, 1, 1, 1;
    1, 1, 1, 1, 0, 0;
    1, 1, 1,1 ,1 ,0;
    1, 1, 1, 1, 1, 1];
%this is not block-lower-triangular
%so matrix elimination must be disabled
config.syn.elimination = 0;

%solve the problem
man = opt_synthesis(sys);
man_channel = opt_synthesis(sys_channel, config);
sol_channel = man_channel.bisect(); 
% 0.9357 with block-triangular
% 0.9341 with relaxed info structure


%% plot solutions
rng(345, 'twister');
M = cell(s-1, 1);
bstar = cell(s-1, 1);
d = 200;

ops_sim = cell(s, 1);
for i = 1:s
    M{i} = rand_quad(d, m(i), L(i));
    bstar{i} = 10*randn(d, 1) + 100*randn(1, 1);
    ops_sim{i} = op_sim_quad(M{i}, bstar{i});
end

tau = 100;
ops_sim{end} = op_sim_l1_hard(tau);


sys_sim = sol_channel.export_sim(ops_sim);

T = 200;
sim = alg_sim(sys_sim, d);
out = sim.sim(T);

plt = alg_plotter(out);
plt.plot({'f', 'res_w', 'res_z'}, 3)

