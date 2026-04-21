function [M_outer] = outer_blkdiag(M1, M2, n1, m1, n2, m2)
%UNTITLED block diagonalization operator for iqc
%utility function for the loop transformation and running cost, allowing
%for flexibility in the lmim routines

M_diag = blkdiag(M1, M2);
P_loop = eye(n1+m1+n2+m2);

ind_perm = [1:n1, (n1+m1) + (1:n2), ...
    n1 + (1:m1), (n1+m1+n2) + (1:m2)];
P_loop(:, ind_perm) = P_loop;
% P_loop = eye(Y

% [E1, E2] = outer_diag_indexer(n1, m1, n2, m2);
M_outer = P_loop* M_diag * P_loop';

end