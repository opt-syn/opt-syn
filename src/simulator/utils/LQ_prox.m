function prox_x = LQ_prox(c, x, n, Q_all, b_all)

    prox_x = zeros(size(x));

    N = length(n);

    for i = 1:N
        ind_curr = (1:n(i)) + sum(n(1:i-1));
        ind_other = setdiff(1:(sum(n)), ind_curr);
    
        %increase monotonicity by adding a local cost function
        % Q_curr = Game{i}.Q;
        % b_curr = Game{i}.b;
        b_self = b_all(ind_curr);

        Q_self = Q_all(ind_curr, ind_curr);
        Q_other= Q_all(ind_curr, ind_other);

        v_self = x(ind_curr);
        v_other= x(ind_other);

        prox_x_curr = LQ_prox_agent(c, v_self, v_other, Q_self, Q_other, b_self);

        prox_x(ind_curr) = prox_x_curr;
    end

end
