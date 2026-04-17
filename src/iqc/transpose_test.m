% z = tf('z', 1);

%how can we do inversion of a system without an inverse?
%warp via bilinear transformation, approximations?

n = 4;
[Af, Bf] = block_fir(n);
Cf = rand(1, n);
Df = 1;
G = ss(Af, Bf, Cf, Df, 1);
Gt = inv(G); 

Gc = d2c(G, 'tustin');
% Gcc = inv(Gc);
Dc = inv(Gc.D);
Gcc = ss(Gc.A - Gc.B*(Dc * Gc.C), Gc.B*Dc, -Dc*Gc.C, Dc);

Gtrec = c2d(Gcc, 1, 'tustin');
% G = (z^3 - z^2 + 2*z - 3)/z^3;
% Ginv = inv(G)