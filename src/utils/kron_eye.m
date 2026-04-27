function [M_lift] = kron_eye(M, c)
%LMIM_KRON_EYE kronecker product with eye(c) when M is an lmim

    if isa(M, 'lmim')
        M_lift = drep( M, c);
        [n, m] = dim(M);
        ind_perm_left = reshape(reshape(1:(n*c), n, c)', [], 1);
        
        P_left = eye(n*c);
        P_left = P_left(:, ind_perm_left);

        ind_perm_right = reshape(reshape(1:(m*c), m, c)', [], 1);
        
        P_right = eye(m*c);
        P_right = P_right(:, ind_perm_right);
        
        M_lift = P_left' * M_lift * P_right;
    else
        M_lift = kron(M, eye(c));
    end

end

