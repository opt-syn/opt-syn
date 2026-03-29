function [p] = model_channel_delay(d1, d2, s)
%MODEL_CHANNEL_DELAY per-channel delays before and after oracle
%evaluation
%
%inputs:
%
%   d1: time delays before the oracle
%   d2: time delays after the oracle
%   s:  number of oracles (if d1 and d2 are scalar)
%
%
%output:
%   p: structure with model (P) and indexers (iz, iy, iu, id

%process 
if nargin == 2
    s = length(d1);
else
    if length(d1) ~= s
        d1 = kron(d1, ones(s, 1));
    end
    if length(d2) ~= s
        d2 = kron(d2, ones(s, 1));
    end
end

%get the pass-through 
p = model_pass_through(s);

G0 = p.P;


    %develop per-channel delays (no crosstalk)

    G1 = ss(zeros(s)); %from u to z
    G2 = ss(zeros(s)); %from w to y   
    
    z = tf('z', 1);
    for i = 1:s
        G1(i, i) =  z^(-d1(i));
        G2(i, i) =  z^(-d2(i));
    end

    %get the final model
    Gd =  G0 * blkdiag(G2, G1);

    for i = 1:length(p.P.A)
        Gd.StateName{i} = sprintf('x%d', i);
    end
    
    %rename the attributes
    for i = 1:(2*s)
        if i <= s
            Gd.InputName{i} = sprintf('w%d', i);
            Gd.OutputName{i} = sprintf('z%d', i);
        else
            Gd.InputName{i} = sprintf('u%d', i-s);
            Gd.OutputName{i} = sprintf('y%d', i-(s));
        end
    end

    %export the model
    p.P = Gd;
    p.opts = alg_options;
end

