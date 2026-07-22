function [Plist] = finite_sum_architecture(s)
%FINITE_SUM_ARCHITECTURE
% finite sum optimization with s functions
%
%this is actually a step towrads an internal model
%
%

Plist = cell(s, 1);

A_base = eye(s);
B_base = zeros(s, 2*s);
C_base = [eye(s); zeros(s)];
D_base = [zeros(s), eye(s); zeros(s), zeros(s)];

n = struct('nz', s, 'nu', s, 'nw', s, 'ny', s, 's', s);

for i = 1:s
    A_curr = A_base;
    A_curr(i, i) = 0;
    B_curr = B_base;
    B_curr(i, i) = 1;
    C_curr = C_base;
    C_curr(i, i) = 0;
    D_curr = D_base;
    D_curr(i+s, i) = 1;

    P_curr = ss(A_curr, B_curr, C_curr, D_curr, 1);
    Plist{i} = genplant(P_curr, n);

end

end

