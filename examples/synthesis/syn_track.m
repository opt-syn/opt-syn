%multi-step method testing
% s = 6;
s = 3;
m = ones(s);
L = (s:-1:1);


theta = pi/8;
Sbeta = blkdiag(1, givens(cos(theta), sin(theta)));
Rbeta = [1, 1, 0];


op_list = cell(s, 1);

for i = 1:s
    op_list{i}= op_sml(m(i), L(i));   
end

sys = opt_system(op_list);
sys.tracking = struct('Sbeta', Sbeta, 'Rbeta', Rbeta);

% config = opt_config();
% config.syn.D_mask = zeros(s);
% config.syn.reduced_order = false;
% man_f = opt_synthesis(sys, config);
% sol_f = man_f.bisect();
% sol_f.rho

% config.syn.D_mask = tril(ones(s), -2);
config = opt_config();
config.syn.D_mask = zeros(s);
config.syn.reduced_order = true;
config.gen.same_rho = true;
config.syn.elimination = false;
% config.gen.same_rho = false;
% config.syn.D_mask(end, 1) = 1;
man = opt_synthesis(sys, config);

sol = man.bisect();
sol.rho