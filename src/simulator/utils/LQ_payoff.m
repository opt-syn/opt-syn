function payoff_opt = LQ_payoff(x_curr, Game)
    
  
    N = length(Game);
    payoff_opt = zeros(N, 1);
    for i = 1:N
        Q_curr = Game{i}.Q;
        b_curr = Game{i}.b;
        c_curr = Game{i}.c;

            payoff_curr = 0.5*(x_curr'*Q_curr*x_curr) + b_curr'*x_curr + c_curr;
            payoff_opt(i) = payoff_curr;

    end
end