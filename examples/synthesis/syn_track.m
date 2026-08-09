%multi-step method testing

s = 6;
m = ones(s);
L = (1:s)*2;


theta = pi/8;
Sbeta = blkdiag(1, givens(cos(theta), sin(theta)));
Rbeta = [1, 1, 0];


op_list = cell(s, 1);

for i = 1:s
    op_list{i}= op_sml(m(i), L(i));   
end


sys = opt_system(op_list);
sys.tracking = struct('Sbeta', Sbeta, 'Rbeta', Rbeta);


config = opt_config();
% config.syn.D_mask = tril(ones(s), -2);
config.syn.D_mask = zeros(s);
% config.syn.D_mask(end, 1) = 1;
man = opt_synthesis(sys, config);

sol = man.bisect();
sol.rho