function [DB] = dhd_basis_0(n)
%DHD_BASIS basis matrices for the set of DHD operators
%
DB_side = dhd_basis(n+1);

DB = DB_side(1:n, 1:n, :);



end

