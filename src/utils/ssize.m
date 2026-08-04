function [sz1, sz2] = ssize(M, ind)
%SSIZE size for lmim types is dim. this is a hack.

if isnumeric(M)
    f = @size;
else
    f = @dim;
end

if nargin == 1
    [sz1, sz2] = f(M);
else
    sz1 = f(M, ind);
    sz2 = [];
end


end

