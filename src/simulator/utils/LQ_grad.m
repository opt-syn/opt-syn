function grad_x = LQ_grad(x, Q_all, b_all)

    grad_x = Q_all*x + b_all;

end
