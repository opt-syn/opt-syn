%synthesis of fixed point algorithms 
%
%set-valued maps (nash equilibria)
    

%% generate parameters of the game
rng(40, 'twister')

c=1;

CONS = true;

if CONS
    s=2;
else
    s = 1;
end

mu = 1.4785;
beta = 0.1605;

BOX = 5;

op1 = op_gen();
op1.monotone = mu;
op1.cocoercive = beta;

if ~CONS
    ops = {op1};
else
    op2 = op_sml(0, inf, 1);
    
    ops = {op1, op2};
end

 

% d1 = [1, 0];
% d2 = [1, 0];
d1 = 0;
d2 = 0;
network = bridge_channel_delay(d1, d2, c);

sys = opt_system(ops,  [], []);


% rho = 0.95;
% rho = 0.7;
% rho = 0.9;
% rho = 0.85;
% rho = 1;
rho = 0.995;
% rho = 0.99
perf_stab = spec_stability(rho);
specs = {perf_stab};

config = opt_config();
config.syn.elimination = true;
config.syn.reduced_order = true;
config.tol.GX_max = 100;
config.tol.GY_max = 100;
if CONS
    config.syn.D_mask = [1, 0; 1, 1]; 
    % config.syn.D_mask = [1 0; 1, 1];
    
else
    config.syn.D_mask = 0; %infeasible
    % config.syn.D_mask = 1; %feasible
end
config.recovery.blocks = true;
man = opt_synthesis(sys, config);

BISECT = true;
    b_opts = bisect_opts;

order = {1, 1};

SOLVE_TYPE = 1;
iqc_warm = sol_ana.cert.iqc_op; 
if SOLVE_TYPE == 2
    b_opts = bisect_opts;
    Niter = 2;
    b_opts.backoff = 1e-3;
    order = {2, [2, 2]};

    [sol_history, vr_history, success] = man.alternate(Niter, order, [], specs, b_opts);

    sol = sol_history{1, end};

elseif SOLVE_TYPE == 1

    [sol_best, v_range] = man.bisect(iqc_warm, specs, b_opts)
    sol = sol_best;
else
    sol = man.solve_single(iqc_warm, specs)
end
