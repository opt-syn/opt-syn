function [sK, sCore] = static_sub(b, s, p)
%STAT Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3 &&  ~isempty(p.opts.tracking)
    Sbeta = p.opts.tracking.Sbeta;
    Rbeta = p.opts.tracking.Rbeta;
else
    Sbeta = 1;
    Rbeta = 1;
end
d = size(Sbeta, 1);


reg = p.regulator;
Gam1 = reg.Gam(:, 1:d);

%in this experiment, only deal with Phi1 = 0 and Phi2 = N.

% Gam1 = reg.Gam(:, 1:d);

if isstruct(b) || istable(b)
    b0 = table2array(b(:, 1:d))';
    b1 = b.b1;
    b2 = b.b2;
else
    b0 = b(1:d);
    b1 = b(end-1);
    b2 = b(end);
end



bz = zeros(1, s);
bz(1) = 1;
bz2 = zeros(1, s);
bz2(end) = 1;
D = [b0 * ones(1, s);  ones(s, 1)*bz*b1 + bz2'*ones(1, s)*b2];

if isa(b, 'genmat')
    sCore = genss(D);
else
    sCore = ss(D);
end
sCore.Ts = 1;

sTrack = ss(Sbeta, [eye(d, d), zeros(d, s)], -Gam1, [zeros(s, d), eye(s)], 1);

sK = sTrack * sCore;
 

end

