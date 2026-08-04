function so = drep_scalar(s, d)
%DREP_SCALAR Summary of this function goes here
%   Detailed explanation goes here

so = zeros(d);
for i = 1:d
    so(i, i) = so(i, i) + s;
end

end

