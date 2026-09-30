%well-posedness check for douglas rachford%analysis of douglas-rachford algorithm
m1 = 1;
b1 = 0.1;

m2= 2;
L2 = 10;

a = sdpvar(4, 1);

P11 = [-2*m1, 1;
     1, 0];
P12 = [0, 1;
       1, -2*b1];

J = [0, 1; 1, 0];
Psi = [L2/(L2-m2), -1/(L2 - m2);
     -m2, 1];
P2 = Psi' * J * Psi;

P3 = J;

Q = blkdiag(P11(1, 1) * a(1) + P12(1, 1) * a(2), P2(1, 1) * a(3), P3(1, 1) * a(4));
S = blkdiag(P11(1, 2) * a(1) + P12(1, 2) * a(2), P2(1, 2) * a(3), P3(1, 2) * a(4));
R = blkdiag(P11(2, 2) * a(1) + P12(2, 2) * a(2), P2(2, 2) * a(3), P3(2, 2) * a(4));


P = [Q, S; S', R];

%douglas-rachford
gamma = 0.3;    %stepsizes
lambda = 0.2;
sK = ss([1], [-gamma*lambda, -gamma*lambda, -gamma*lambda], [1; 1; 1],...
    [-gamma, 0, 0; -gamma, 0, 0; -2*gamma, -gamma, -gamma],1);

E = sK.D;


right = [E; eye(3)];

Pi = right' * P * right;

cons = [a >= 0, sum(a)==1, Pi <= -eye(3)*1e-3];

opts = sdpsettings('solver', 'sdpt3');
optimize(cons, norm(a, 2),  opts)
