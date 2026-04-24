function [Plist] = finite_sum_primitives(s)
%FINITE_SUM_PRIMITIVES
% finite sum optimization with s functions
%
%this is actually an internal model
%
%

Plist = cell(s, 1);

A_base = [];
B_base = zeros(0, s+1);
C_base = zeros(s+1, 0);
D_base = zeros(s+1);

n = struct('nz', s, 'nu', 1, 'nw', s, 'ny', 1, 's', s);

for i = 1:s
    A_curr = A_base;
    B_curr = B_base;
    C_curr = C_base;
    D_curr = D_base;
    D_base(end, i) = 1;
    D_base(i, end) = 1;

    P_curr = ss(A_curr, B_curr, C_curr, D_curr, 1);
    Plist{i} = genplant(P_curr, n);

end

end

