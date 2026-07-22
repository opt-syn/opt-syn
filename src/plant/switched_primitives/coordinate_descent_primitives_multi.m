function [Plist, bin_seq] = coordinate_descent_primitives_multi(c, simul)
%COORDINATE_DESCENT_PRIMITIVES_MULTI
% block coordinate descent with c blocks
% where 'simul' entries are updated at a time (all combinations)

if nargin == 1
    simul = 1;
end


simul = reshape(simul, 1, []);

bin_seq = (dec2bin(0:2^c-1)' - '0')';
bin_sum = sum(bin_seq, 2);

ind_active = any(bin_sum == simul, 2);

bin_seq = bin_seq(ind_active, :);

nbin = sum(ind_active);

Plist = cell(nbin, 1);


A_base = eye(c);
% A(1, 1) = 0;
% B = [[0, 1]; zeros(c-1, 2*c)];
B_base = zeros(c, 2*c);
% B(1, c+1) = 1;
% C = [A; zeros(c)];
C_base = [eye(c); zeros(c)];
D_base = [zeros(c), zeros(c); eye(c), zeros(c)];
% D(1, c+1) = 1;

% delay_curr = ss(A, B, C, D, 1);

n = struct('nz', c, 'nu', c, 'nw', c, 'ny', c, 's', 1);






% ind_odd = 2*(0:c-1)+1;
% ind_even = ind_odd+1;
for j = 1:nbin
    bin_curr = bin_seq(j, :);
    A_curr = A_base;
    B_curr = B_base;
    C_curr = C_base;
    D_curr = D_base;
    for i = 1:c
        if bin_curr(i)
            A_curr(i, i) = 0;    
            B_curr(i, i+c) = 1;    
            C_curr(i, i) = 0;    
            D_curr(i, i+c) = 1;
        
        end
    end

    P_curr = ss(A_curr, B_curr, C_curr, D_curr, 1);
    Plist{j} = genplant(P_curr, n);

end
end

