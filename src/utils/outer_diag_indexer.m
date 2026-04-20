function [E1, E2] = outer_diag_indexer(n1, m1, n2, m2)
%UNTITLED block diagonalization operator for iqc
%utility function for the loop transformation and running cost, allowing
%for flexibility in the lmim routines

E1 = [eye(n1, n1), zeros(n1, m1);
    zeros(n2, m1), zeros(n2, m1);
    zeros(n1, m1), eye(n1, m2);
    zeros(n2, m1), zeros(n2, m2)];

E2 = [zeros(n1, n1), zeros(n1, m1);
    eye(n2, n2), zeros(n2, n1);
    zeros(m2, m2), zeros(m2, m1);
    zeros(m1, m2), eye(m1, m1)];


end