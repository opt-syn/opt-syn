function [T, Ti, G, Gi] = rec_similarity(vars, p)
%REC_similarity Summary of this function goes here
%   Detailed explanation goes here


%need to adjust this for multi-input multi-output cases
iz = p.iz;
iw = p.iw;
iy = p.iy;
iu = p.iu;

%recover the transformation matrix

if isfield(vars, 'S')
    n = length(vars.S);
    S = (vars.S);
    Y = vars.Y;
    X = vars.X;
else
    n = length(vars.P);
    X = vars.P;
    Y = vars.H;
    S = eye(n);
end

%get left-side entries


I = eye(n);
J = S - Y * X;
[Up, Sig, Vp] = svd(J);

% U = Up*Sig;
ssig = sqrt(Sig);
srsig = diag(1./(diag(ssig)));

U = ssig*Vp';
V = Up*ssig;

Uinv = Vp*srsig;
Vinv = srsig*Up';


%get right-side entries

Z1 = (Vinv*(I - X * Y')')';
Z2 = (Vinv* (-U * Y')')';

Z34 = [X, Z1; U, Z2] \ [zeros(n); eye(n)];

Z3 = Z34(1:n, :);
Z4 = Z34((n+1):end, :);

T = [eye(n), Y'; zeros(n), V'];
Ti = [eye(n), -Y' * Vinv'; zeros(n), Vinv'];

G = [X, Z1; U, Z2];
Gi = [Y', Z3; V', Z4];

end

