NL = 3;
Llist = logspace(-2, 3, NL);

orderlist = {{[3, 1], [3, 1]}
    {[3, 0], [3, 0]}};
Norder = length(orderlist);

rho = zeros(NL, 2);



%analysis of proximal heavy ball algorithm




for i = 1:NL
    m= 1;
    L = Llist(i);

    %different operator classes for op1
    op1_sml = op_sml(m, L);
    op1_quad = op_quad(m, L);
    
    %indicator function for op1
    op2 = op_sml(0, inf, 1);
    %proximal heavy ball
    alpha = 2/(m+L);
    beta = 0.65;
    
    A = [1+beta, -beta; 1, 0];
    B = [-alpha, -alpha; 0 , 0];
    C = [1, 0; 1+beta, -beta];
    D = [0, 0; -alpha, -alpha];
    
    sK = ss(A, B, C, D, 1);

    sys_sml = opt_system({op1_sml,  op2}, [], sK);
    sys_quad = opt_system({op1_quad, op2}, [], sK);

    man_sml = opt_analysis(sys_sml);
    man_quad = opt_analysis(sys_quad);

    for j = 1:orderlist
        sol_sml = man_sml.bisect(orderlist{j});
        rho_sml(i, j) = sol_sml.rho;
    
        sol_quad = man_quad.bisect(orderlist{j});
        rho_quad(i, j) = sol_quad.rho;
    end

end

% [0.8971, sol_quad.rho]
% [sol_sml.rho, sol_quad.rho]