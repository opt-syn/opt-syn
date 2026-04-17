function [M_lift] = kron_eye(M, c)
%LMIM_KRON_EYE kronecker product with eye(c) when M is an lmim

    if isa(M, 'lmim')
        M_lift = drep( M, c);
        [n, m] = dim(M);
        ind_perm = reshape(reshape(1:(n*c), n, c)', [], 1);
        
        P = eye(n*c);
        P = P(:, ind_perm);
        
        M_lift = P' * M_lift * P;
    else
        M_lift = kron(M, eye(c));
    end

end

