function [xproj] = weighted_l1_proj(y, rad, w)
%WEIGHTED_L1_PROJ projection onto weighted l1 ball
%
%used for prox operators with inhomogenously-weighted lasso
if nargin < 3
    w = ones(size(x));
end

s = sign(y);

zu = abs(y)./w;

[z, perm] = sort(zu);

wy = cumsum(w(perm).*abs(y(perm)), 'reverse');
w2 = cumsum(w(perm).^2, 1, 'reverse');

%algorithm 2 https://www.sciencedirect.com/science/article/pii/S0004370222000236

th1 = ((wy - rad) ./ w2) - z;

j = find(th1 > 0, 1, 'last');

if isempty(j)
    j = 0;
end
lam = (wy(j+1) - rad)/w2(j+1);


xproj = sign(y) .* max(abs(y) - w* lam, zeros(size(y)));

end

