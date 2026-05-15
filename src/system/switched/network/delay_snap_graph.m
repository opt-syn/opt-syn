function [Gcon] = delay_snap_graph(N1, N2)
%DELAY_SNAP_GRAPH snapping logic (packet drop)
%   Detailed explanation goes here

if nargin < 5
    toggle = 0;
end

% A1 = zeros(N1+1);
% A1(1, :) = 1
A1 = [ones(N1, 1), eye(N1); 1, zeros(1, N1)];
A2 = [ones(N2, 1), eye(N2); 1, zeros(1, N2)];
% A2 = zeros(N2 + 1);

Gcon = kron(A1, A2);



%TODO: fix the snapping logic
% if snap
%     Gcon(end, end) = 0;
% end



