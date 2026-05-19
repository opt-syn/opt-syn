function [sys_warp, M, T] = warp_system(sys, alpha, beta)
%WARP_SYSTEM use a mobius transformation to warp a system via 
% z -> (alpha z - b)/(alpha z + b)


if nargin == 1
    alpha = 2;
    beta = -2;
end

[A1, B1, C1, D1] = ssdata(sys);


[n1, m1] = size(B1);
p1= size(C1, 1);
M = kron([alpha, -beta; alpha, beta], eye(n1));

T = [alpha*A1 + beta*eye(n1), alpha*B1;
    zeros(m1, n1), eye(m1)];

ABM = [A1, B1; eye(n1), zeros(n1, m1)];

ABM_warp = M * ABM * inv(T);

CDM = [C1, D1;
    zeros(p1, n1), eye(p1)] * inv(T);



CDM_warp = CDM * inv(T);

Aw1 = ABM_warp(1:n1, 1:n1);
Bw1 = ABM_warp(1:n1, (n1+1):end);

Cw1 = CDM_warp(1:p1, 1:n1);
Dw1 = CDM_warp(1:p1, (n1+1):end);

sys_warp = ss(Aw1, Bw1, Cw1, Dw1, 0);


end