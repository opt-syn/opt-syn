%generate the state space subsystems
P =cell(4, 1);
nx = 2;
P{1} = [0 , 0 , 1 , 0;
    0 , 0.2 , 0 , 1; 
    0 , 1 , 0 , 0 ; 
    -1 , 0 , 0 , 0];

P{2} = [0.2 , 0 , 0.25 , 0;
    0 , 0.9 , 0 , 1 ; 
    0 , 1 , 0 , 0 ; 
    -0.4 , 0 , -0.5 , 3];

P{3} = [-0.3 , 0 , 0.5 , 0;
    0 , -0.5 , 0 , 1 ;
    0 , 1 , 0 , 0 ; 
    -0.3 , 0 , 0.5 , 1];

P{4} = [1.2 , 0 , 1 , 0;
    0 , -0.2 , 0 , 1 ;
    0 , 1 , 0 , 0 ;
    1 , 0 , 0 , 2];

%break up the subsystems into genplants
Plist = cell(4, 1);
n = struct('nw', 1, 'nz', 1, 'ny', 1, 'nu', 1);
for i = 1:4
    A = P{i}(1:nx, 1:nx);
    B = P{i}(1:nx, (nx+1) : end);
    C = P{i}((nx+1) : end, 1:nx);
    D = P{i}((nx+1) : end, (nx+1) : end);
    Plist{i} = genplant(ss(A, B, C, D, 1), n);
end

network = genplant_poly(Plist);

%define the operator
m = 1; L = 2;
ops = {op_sml(m, L)};

sys = opt_system_periodic(ops, network, []);

%only allow gradients
config =opt_config();
config.syn.prox = 0;

%pose and solve
man= opt_synthesis(sys, config);
sol= man.bisect();   %0.7028

%% begin simulation
d = 50;
Q = rand_quad(d, m, L);
bstar = randi(101, [d, 1]) - 50;
op1 = op_sim_quad(Q, bstar);

ops_sim = {op1};

sys = sol.sys.export_sim(ops_sim);

T = 60;

sim  = alg_sim(sys, d);
nx = sys.get_alg(1).nx;
sim_out = sim.sim(T);

plt = alg_plotter(sim_out);
plt.plot({'x', 'w', 'res_w', 'mode', 'z', 'f'}, 1)


%% regulator equations
sim_out_long = sim.sim(6*T);
dstar = -sim_out_long.z(:, :, end);
plt = plt.add_opt_sig(sol.regcl, dstar);

%% new
plt = plt.add_opt_sig(sol.regcl, dstar);
plt.plot_4_err(2);    %plot the tracking error
plt.plot_4_sq_err(3); %plot the  squared norm of the tracking error
