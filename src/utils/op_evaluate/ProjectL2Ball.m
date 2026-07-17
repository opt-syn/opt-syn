function [xproj] = ProjectL2Ball(x, ballRadius, x0)
%PROJECTL2BALL Projection onto the L2 ball (unit sphere)

if nargin < 3
    x0 = 0;
end

if norm(x-x0, 2) < ballRadius
    xproj = x;
else
    dir = (x-x0)/norm(x-x0, 2);

    xproj = x0 + dir*ballRadius;
end

end

