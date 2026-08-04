function [X] = ball_sample(N, d)
%BALL_SAMPLE Randomly sample N points from a unit d-dimensional ball
%d=3 sphere is a 2-sphere. S^(d-1)
%
%Input:
%   N:  Number of points to sample
%   d:  Dimension of ball
%
%Output:
%   X:  Points on ball

%dropped coordinate method
u = sphere_sample(N, d+2);
X = u(:, 1:d);

end

