% syms L m rho G
% syms G
% assume([L, m, rho, G], 'real')
CONS = true;
% STAB = true;
STAB = false;
G = lmim('G', 1, 1);
% G =  3.0262;

L = 10;
m = 1;
% rho = 0.92;
rho = 0.85;
rhoi = 1/rho;
gam = 2/(L+m);
% gam = 1/L;

sig = 1/(L-m);
k = (L-m)/(L+m);

A = rhoi*k;
B = -rhoi*gam*[1, 1];
C = [1; k];
D = [-sig, 0;
    -gam, -gam];

P = ss(A, B, C, D, 1);
 

AB = [1, 0, 0;
 A, B];

IA = [1; A];
Gblock = [G 0 ; 0 -G];
Ablock = AB' * Gblock * AB;

CD = [zeros(2, 1), eye(2)
      C, D];

%strict negative realness
% epass = 1e-4;
% epass = 0;
epass = 1e-3;
Cblock = -CD' * kron([0, 1; 1, epass], eye(2)) * CD;



if STAB
    K = IA'* Gblock * IA;
else
    K = Ablock + Cblock;
end

cons = [];
cons = append_lmi(cons, G - 1e-4, 1);
cons = append_lmi(cons, K, 1);

if CONS
    s = trace(G);
    cons = lmis(cons, s, 'c');   
end
[lmi_out,info_out]=lmisolve(cons);

if ~lmi_out.status

    G_rec = double(double(G, lmi_out))
    K_rec = double(double(K, lmi_out))
    % eK = eig(K_rec)'
else
    disp('infeasible')
end