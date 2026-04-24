function [p] = model_general(nconn, s, eigset)
%MODEL_GENERAL create a random model for the network dynamics
%   feedback and crosstalk is allowed
%
%
%Inputs:
%   nconn:  number of states in the model
%   s:      number of oracles
%   eigset: assign the maximal eigenvalue of the system

p = struct;

p.iz = 1:s;
p.iy = s + (1:s);
p.id = 1:s;
p.iu = s + (1:s);

Gd = drss(nconn, 2*s, 2*s);
Gd.Ts = 1;   

if nargin == 3 && nconn>0
    Gd.A = eigset* Gd.A / max(abs(eig(Gd.A)));
end
    
    
for i = 1:(2*s)
    if i <= s
        Gd.InputName{i} = sprintf('w%d', i);
        Gd.OutputName{i} = sprintf('z%d', i);
    else
        Gd.InputName{i} = sprintf('u%d', i-s);
        Gd.OutputName{i} = sprintf('y%d', i-(s));
    end
end

for i = 1:nconn
    Gd.StateName{i} = sprintf('x%d', i);
end

  %try to promote a solution to the regulator equation 
    %because random generation might not work
    ind_miss = sum([Gd.B(:, p.iu); Gd.D(p.iz, p.iu)]~=0, 1)==0;
    for i = 1:s
        if ind_miss(i)
            Gd.D(i, i+s) = 1;
        end
    end

p.P = Gd;


p.iz = 1:s;
p.iy = s + (1:s);
p.iw = 1:s;
p.iu = s + (1:s);

p.opts = alg_options;

end

