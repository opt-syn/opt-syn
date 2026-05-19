function [sys_warp] = warp_system_mat(sys, alpha, beta)
%WARP_SYSTEM use a mobius transformation to warp a system via 
% z -> (alpha z - b)/(alpha z + b)


if nargin == 1
    alpha = 2;
    beta = -2;
end

[A, B, C, D] = ssdata(sys);
n = size(A);
ad = alpha*A - beta*eye(n);
ap = inv(alpha*A + beta*eye(n));



Aw = ad * ap;
Bw = 2*alpha*beta*ap*B;
Cw = C*ap;
Dw = D - alpha*Cw*B;

sys_warp = ss(Aw, Bw, Cw, Dw, 0);

end