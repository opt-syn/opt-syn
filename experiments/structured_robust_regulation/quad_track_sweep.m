%tracking a quadratic with a rotating optimal solution

%define the oracles
s = 4;
m = [1, 0, 2, -2];
L = [3, 10, 3, -1];

ops_sim = cell(s, 1);

for i = 1:s
    ops_sim{i}= op_sml(m(i), L(i));   
end


%define the rotation of the optimal solution

theta = pi/8;
Sbeta = blkdiag(1, givens(cos(theta), sin(theta)));
Rbeta = [1, 1, 0];

%form the specification (noisy gradient measurements)
network = bridge_pass_through(s);

network = network.add_oracle_input(1:s, []);
% network = network.perf_output_z(1:s);
network = network.perf_output_con;

GAIN = 100;
spec_gain = spec_e2e(GAIN, 1:s, 1:s);
spec_gain.target = true;
rho = 0.95;

spec_stab = spec_stability(rho);
specs = {spec_stab, spec_gain};

%form the system

sys = opt_system(ops_sim, network);
tracking = struct('Sbeta', Sbeta, 'Rbeta', Rbeta);
sys.tracking = tracking;

%solve the problem (bisection)
config = opt_config();
% config.syn.D_mask = zeros(s);
config.syn.D_mask = tril(ones(s), -1);
config.gen.same_rho = true;

man = opt_synthesis(sys, config);

sol = man.solve_single([], specs);


%% simulate and plot
rng(32, 'twister');

d = 5;

%tracking of the solution
SA = kron(Sbeta, eye(d));
SAy = kron(Rbeta, eye(d));
eta0 = [randi(201, [d, 1]) - 100; 
    (randi(31, [2*d, 1]) + 60).*sign(2*rand(2*d, 1) - 1)];

shift = @(k) SAy * (SA^k) * eta0;
shift0 = @(k) 0*SAy * (SA^k) * eta0;

M = cell(s, 1);
bstar_center = cell(s, 1);
bstar = cell(s, 1);
bstar_0 = cell(s, 1);
ops_sim = cell(s, 1);
ops_sim_0 = cell(s, 1);

for i = 1:s
    M{i} = rand_quad(d, m(i), L(i));
    bstar_center{i} = randi(21, [d, 1]) - 10;    
    bstar{i} = @(k) bstar_center{i} - shift(k);
    bstar_0{i} = @(k) - shift(k);
    
    ops_sim{i} = op_sim_quad(M{i}, bstar{i});  
    ops_sim_0{i} = op_sim_quad(M{i}, bstar_0{i});  
end

%simulate

sys_sim = sol.export_sim(ops_sim);
sim = alg_sim(sys_sim, d);
sim.sampler.x0 = 7*(2*rand(sys_sim.n, d)-1);
wp_std = 10;
wp_alpha = -0.1;
% sim.sampler.wp = @(k, param) wp_std*rand(1, s*d);
sim.sampler.wp = @(k, param) wp_std*(2*rand(s, d) - 1) * exp(k*wp_alpha);


T = 100;
sim_out= sim.sim(T);
plt = alg_plotter(sim_out);

plt.plot({'w', 'wp' 'z', 'zp'}, 1)


%check the gain
wp_sq = squeeze(sum(sim_out.wp.^2, [1, 2]));
zp_sq = squeeze(sum(sim_out.zp.^2, [1, 2]));
rhod = rho.^(-2*(0:T-1))';
gain = cumsum(rhod .* zp_sq) ./ cumsum(rhod .* wp_sq);


%% connect to an external system

gain_alg = sqrt(sol.objective);

z = tf('z', 1);
Pext_1 = minreal(ss(0.3*z^(-1) - 0.2*z^(-2) -0.1 * z^(-3)));
% Pext_1 = 0.5*minreal(ss(0.5 - 0.5*z^(-1) - 0.3*z^(-2) -0.2 * z^(-3)));

Pext = ss_kron_eye(Pext_1, s);
% Pext = 0.5*eye(s);

Prext = rhotrafo(Pext, rho);
hrext = hinfnorm(Prext);

small_gain = hrext * gain_alg


Pall = lft(lft(sys_sim.P.P, sys_sim.K.P), Pext);

sys_ext = opt_system(ops_sim, [], Pall);
sys_ext.tracking = tracking;


sim_ext = alg_sim(sys_ext, d);
% sim_ext.sampler.x0 = 7*(2*rand(sys_ext.n, d)-1);

T = 300;
sim_out_ext= sim_ext.sim(T);
plt_ext = alg_plotter(sim_out_ext);
plt_ext.plot_6f(2);


reg_ext = regulator_lti(sys_ext);
regcl = reg_ext.check_regulator();


%% get properties of the external system
Nrho = 100;
rhol = linspace(0.8, 1, Nrho);
gain_ext = zeros(Nrho, 1);
for i = 1:Nrho
    gain_ext(i) = hinfnorm(rhotrafo(Pext, rhol(i)));
end
figure(4)
clf
plot(rhol, gain_ext)
