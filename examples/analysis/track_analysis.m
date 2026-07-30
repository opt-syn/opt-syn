%multi-step method testing


m = [1, 1, 1, 1];
L = [2, 4, 6, 8];
s = length(m);

theta = pi/8;
Sbeta = blkdiag(1, givens(cos(theta), sin(theta)));
Rbeta = [1, 1, 0];


op_list = cell(s, 1);

for i = 1:s
    op_list{i}= op_sml(m(i), L(i));   
end


%static parameters
% b0 = [-0.0600696558727804
%     -0.116093656856806
%     -0.00543012633720185]; 
% 
% b1 = -0.00543012633720185;
% b2 = -0.116093656856806;

b0 = [-0.05; -0.1; -0.05];
b1 = -0.05;
b2 = -0.1;


AK = Sbeta;
BK = b0 * ones(1, s);
CK = ones(s, 1) * Rbeta;

DK = zeros(s);
DK(end, :)= DK(end, :) + b1;
DK(:, 1)= DK(:, 1) + b2;


K = ss(AK, BK, CK, DK, 1);


sys = opt_system(op_list, [], K);


man = opt_analysis(sys);

order = {[1, 0], [1, 0], [1, 0], [1, 0]};
sol = man.bisect(order);
sol.rho