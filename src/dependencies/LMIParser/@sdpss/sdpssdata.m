function [A, B, C, D] = sdpssdata(G)
%SDPSSDATA Summary of this function goes here
%   Detailed explanation goes here

if isa(G, 'sdpss')
    [A, B, C, D] = ssdata(G);
else
    A = G.A;
    B = G.B;
    C = G.C;
    D = G.D;
end


end

