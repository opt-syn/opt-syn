function [DB] = dhd_basis(n)
%DHD_BASIS basis matrices for the set of DHD operators
%
N = n*(n-1)/2;
prot = [1, -1; -1, 1];

DB = zeros(n, n, N);
count = 1;
for i = 1:n
    for j = 1:(i-1)

        E = zeros(n, 2);
        E(i, 1) = 1;
        E(j, 2) = 1;

        DE = E * prot * E';

        DB(:, :, count) = DE;
        
        count = count + 1;
    end
end


end

