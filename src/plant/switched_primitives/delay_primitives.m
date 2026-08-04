function [Plist, Gcon] = delay_primitives(d1, d2, r1, r2)
%DELAY_INT_SYSTEMS Generate all subsystems for time-varying delays
%along with their switching transition matrix
%Input:
%   d1: delay before the oracle
%   d2: delay after the oracle
%   r1: change in delay before the oracle
%   r2: change in delay after the oracle
%
%Output:
%   Plist:  cell array of subsystems (with compatible dimensions)
%   Gcon:   graph representing the connectivity between subsystems

N1 = length(d1);
N2 = length(d2);

Nr1 = length(r1);
Nr2 = length(r2);

Plist = cell(N1, N2);

NN = N1 * N2;

d1max = max(d1);
d2max = max(d2);

%generate the subsystems
for i = 1:N1    
    for j = 1:N2
        d1_curr = d1(i);
        d2_curr = d2(j);

        Pcurr = delay_subsystem(d1max, d2max, d1_curr, d2_curr);
        n = struct('nw', 1, 'nz', 1, 'ny', 1, 'nu', 1);
        
        P = genplant(Pcurr, n);
        Plist{i, j} = P;
    end
end

%generate the switching graph
[Gcon] = delay_switch_graph(N1, N2, r1, r2);
