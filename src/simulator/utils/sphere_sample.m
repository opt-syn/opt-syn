function [X] = sphere_sample(N, d)
%SPHERE_SAMPLE Randomly sample N points from a unit d-dimensional sphere
%d=3 sphere is a 2-sphere. S^(d-1)
%
%Input:
%   N:  Number of points to sample
%   d:  Dimension of sphere
%
%Output:
%   X:  Points on sphere

%dropped coordinate method
U = randn(N, d);
normU = sqrt(sum(U.^2, 2));
X = U ./ normU;

end

