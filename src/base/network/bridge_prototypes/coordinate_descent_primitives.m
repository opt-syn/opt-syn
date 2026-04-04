function [Plist] = coordinate_descent_primitives(c)
%COORDINATE_DESCENT_SYSTEM 
% block coordinate descent with c blocks

Plist = cell(c, 1);

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
for i = 1:c
    A_curr = A_base;
    A_curr(i, i) = 0;
    B_curr = B_base;
    B_curr(i, i+c) = 1;
    C_curr = C_base;
    C_curr(i, i) = 0;
    D_curr = D_base;
    D_curr(i, i+c) = 1;

    P_curr = ss(A_curr, B_curr, C_curr, D_curr, 1);
    Plist{i} = bridge(P_curr, n);

end

end

