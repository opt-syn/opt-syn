function [p] = model_pass_through(s)
%create a network model in which there are no dynamics
%
%Input:
%   s: number of oracles
%   p: structure with model (P) and indexers (iz, iy, iu, id
p = struct;

p.iz = 1:s;
p.iy = s + (1:s);
p.iw = 1:s;
p.iu = s + (1:s);

G0 = ss([zeros(s), eye(s); eye(s), zeros(s)]);
    
    
for i = 1:(2*s)
    if i <= s
        G0.InputName{i} = sprintf('w%d', i);
        G0.OutputName{i} = sprintf('z%d', i);
    else
        G0.InputName{i} = sprintf('u%d', i-s);
        G0.OutputName{i} = sprintf('y%d', i-(s));
    end
end

p.P = G0;
p.opts = alg_options;

end