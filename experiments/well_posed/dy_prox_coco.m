%well-posedness check for douglas rachford%analysis of douglas-rachford algorithm
m1 = 1;
b1 = 0.1;

a = sdpvar(2, 1);

P11 = [-2*m1, 1;
     1, 0];
P12 = [0, 1;
       1, -2*b1];

Q = P11(1, 1) * a(1) + P12(1, 1) * a(2);
S = P11(1, 2) * a(1) + P12(1, 2) * a(2);
R = P11(2, 2) * a(1) + P12(2, 2) * a(2);


P = [Q, S; S', R];

%douglas-rachford

E = 0;
right = [E; eye(1)];

% Pi = right' * P * right;
% 
% cons = [a >= 0, sum(a) == 1, Pi <= -eye(1)*1e-3];
cons = [a >= 0, sum(a) == 1];
% 
obj=norm(a, 2)^2;
opts = sdpsettings('solver', 'sdpt3');
optimize(cons, obj,  opts)

av = value(a);
