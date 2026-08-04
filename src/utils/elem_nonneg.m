function [cons] = elem_nonneg(M, cons, LMILAB)
%ELEM_NONNEG enforce that all entries of the lmim matrix M are nonnegative



if nargin < 2
    cons= [];
end

if nargin < 3
    LMILAB = 1;
end

if ~isempty(M)
    
%used for the LMILAB interface

if LMILAB
    [n, m] = dim(M);
else
    [n, m] = size(M);
end

    cv = cell(n*m, 1);
    count = 1;
    for i =1:n
        e = zeros(1, n);
        e(i) = 1;
        for j = 1:m
            ej = zeros(m, 1);
            ej(j) = 1;
        
            
            vc = e * M * ej;
            if ~isnumeric(vc)
                cv{count} = vc;
            end
            count = count + 1;
        end
    end

    cons = append_lmi(cons, cv, LMILAB);

end

end