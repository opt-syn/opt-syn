function [Gss] = ss(G)
%SS turn into a state space system

if isa(G, 'sdpss')
    [A, B, C, D] = ssdata(G);
    Gss = sdpss(A, B, C, D, 1);
else
    G
end


end

