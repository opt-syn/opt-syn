function [M] = rand_quad(d, m, L)
%RAND_QUAD a randomly generated symmetric matrix with eigenvalues
%between m and L.
%
%
%Args
%   d:  dimension of the matrix
%   m:  minimum eigenvalue of matrix
%   L:  maximum eigenvalue of matrix
%
%Returns:
%   M:  the matrix


    if m == L
        M = eye(d)* m;
    else
        e = m+[0; rand(d-2, 1); 1]*(L-m);
        [Q, ~] = qr(randn(d));    
        M = Q'*diag(e)*Q;
    end
end

