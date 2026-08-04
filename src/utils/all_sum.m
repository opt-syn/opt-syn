function [MS] = all_sum(M, LMILAB)
%ALL_SUM Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    LMILAB = 1;
end

%used for the LMILAB interface

if LMILAB
    [n, m] = dim(M);
else
    [n, m] = size(M);
end

MS = ones(1, n) * M * ones(m, 1);

end

