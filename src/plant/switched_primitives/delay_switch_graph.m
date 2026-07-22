function [Gcon] = delay_switch_graph(N1, N2, r1, r2)
%DELAY_SWITCH_GRAPH graph of the switching
%   Detailed explanation goes here

snap = 0;
if nargin < 5
    toggle = 0;
end

NN = N1*N2;
Gcon = zeros(NN, NN);

for i = 1:N1
    i_stagger = r1 + i;
    i_stagger(i_stagger <= 1) = 1;
    i_stagger(i_stagger  >= N1) = N1;
    if snap
        i_stagger = [i_stagger, 1];
    end
    i_stagger = unique(i_stagger);  

    
    for j = 1:N2
        sys_curr = sub2ind([N1, N2], i, j);

        
        j_stagger = r2 + j;       
        j_stagger(j_stagger <= 1) = 1;
        j_stagger(j_stagger  >= N2) = N2;
        if snap
            j_stagger = [j_stagger, 1];
        end
        j_stagger = unique(j_stagger);
     

        cross_j = kron(ones(1, length(i_stagger)), j_stagger );
        cross_i = kron(i_stagger, ones(1, length(j_stagger )));

        ind_destination = sub2ind([N1, N2], cross_i, cross_j);

        Gcon(sys_curr, ind_destination) = 1;
    end
end

%TODO: fix the snapping logic
% if snap
%     Gcon(end, end) = 0;
% end



