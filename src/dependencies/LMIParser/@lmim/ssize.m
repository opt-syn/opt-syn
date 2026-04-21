function [sz] = ssize(M, ind)
%SSIZE size for lmim types is dim. this is a hack.

if isnumeric(M)
    f = @size;
else
    f = @dim;
end

if nargin == 1
    sz = f(M)
else
    sz = f(M, ind)
end


end

