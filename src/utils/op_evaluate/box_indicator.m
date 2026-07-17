function [f] = box_indicator(z, BOX)
%BOX_INDICATOR Summary of this function goes here
%   Detailed explanation goes here

check = all([abs(z)<=BOX + 1e-7]);

if check 
    f = 0;
else
    f = Inf;
end

end

