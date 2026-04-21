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

    for i =1:n
        e = zeros(1, n);
        e(i) = 1;
        for j = 1:m
            ej = zeros(m, 1);
            ej(j) = 1;
        
            vc = e * M * ej;
            cons = append_lmi(cons, vc + vc', 1);
        end
    end

end

end