function [cons] = dhd_impose(M, cons, LMILAB)
%DHD_IMPOSE: impose that a symmetric matrix M is doubly hyper dominant
%(DHD)


if nargin < 2
    cons= [];
end

if nargin < 3
    LMILAB = 1;
end

%used for the LMILAB interface

if LMILAB
    [n, ~] = dim(M);
else
    [n, ~] = size(M);
end

    
    
    M1 = M*ones(n);
    for i =1:n
        e = zeros(1, n);
        e(i) = 1;
        
        %main diagonal terms are dominant       
        cons = append_lmi(cons, e*M1, 1);

        %off-diagonal terms are nonnegative
        for j = 1:i-1
            ej = zeros(n, 1);
            ej(j) = 1;
        
            vc = e * M * ej;
            cons = append_lmi(cons, -vc, 1);
        end
    end

end

