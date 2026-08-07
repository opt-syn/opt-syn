function [Plist] = coordinate_descent_limited(c)
%COORDINATE_DESCENT_LIMITED
% block coordinate descent with c blocks
% but information is only transmitted through for the active block

Plist = cell(c, 1);

A_base = eye(c);
% A(1, 1) = 0;
% B = [[0, 1]; zeros(c-1, 2*c)];
B_base = zeros(c, 2*c);
% B(1, c+1) = 1;
% C = [A; zeros(c)];
C_base = [eye(c); zeros(c)];
D_base = [zeros(c), zeros(c); zeros(c), zeros(c)];

n = struct('nz', c, 'nu', c, 'nw', c, 'ny', c, 's', 1);

i = 1;
    A_curr = A_base;
    A_curr(i, i) = 0;
    B_curr = B_base;
    B_curr(i, i+c) = 1;
    C_curr = C_base;
    C_curr(i, i) = 0;
    D_curr = D_base;
    D_curr(i, i+c) = 1;
    D_curr(i+c, i) = 1;

    P_curr = ss(A_curr, B_curr, C_curr, D_curr, 1);

    Plist = genplant(P_curr, n);

end

