function [Q_list, b_list, c_list] = rand_LQ_game(n, mu_boost);
%RAND_LQ_GAME Summary of this function goes here
%   Detailed explanation goes here
N = length(n);
d = sum(n);
if nargin < 2
    mu_boost = 0;
end

Q_list = cell(N, 1);
b_list = cell(N, 1);
c_list = cell(N, 1);

qscale = 1;

for i = 1:N
    Q_curr = -qscale* rand(d, d);    

    Q_curr = Q_curr + Q_curr';
    
    
    ind_curr = (1:n(i)) + sum(n(1:i-1));
    ind_other = setdiff(1:(sum(n)), ind_curr);

    %increase monotonicity by adding a local cost function
    Q_self = Q_curr(ind_curr, ind_curr);
    Q_self = Q_self + mu_boost*eye(n(i));
    Q_curr(ind_curr, ind_curr) = Q_self;

    Q_other = Q_curr(ind_curr, ind_other);

    b_curr = 50*randn(d, 1);

    b_self = b_curr(ind_curr);

    c_curr = randn(1, 1);

    %global parameters
    Game_curr = struct;
    Game_curr.Q = Q_curr;
    Game_curr.b = b_curr;
    Game_curr.c = c_curr;

    Q_list{i} = Q_curr;
    b_list{i} = b_curr;
    c_list{i} = c_curr;

end
