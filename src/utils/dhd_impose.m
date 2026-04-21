function [cons] = dhd_impose(M, cons, LMILAB)
%DHD_IMPOSE: impose that a matrix M is doubly hyper dominant
%(DHD)


if nargin < 2
    cons= [];
end

if nargin < 3
    LMILAB = 1;
end

%used for the LMILAB interface

if LMILAB
    [n, m] = dim(M);
else
    [n, m] = size(M);
end

    
    
    M1_h = M*ones(m, 1);
    M1_v = ones(1, n) * M;

    
    cons = elem_nonneg(M1_h, cons, LMILAB);
    cons = elem_nonneg(M1_v, cons, LMILAB);

    for i =1:n
        %off-diagonal terms are nonnegative
        ei = zeros(1, n);
        ei(i) = 1;
        for j = 1:m
            if i ~= j
                ej = zeros(m, 1);
                ej(j) = 1;
            
                vc = ei * M * ej;
                cons = append_lmi(cons, -vc, 1);
            end
        end
    end

end

