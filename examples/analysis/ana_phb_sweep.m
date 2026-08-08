%parameters of the sweep
NL = 200;
Llist = logspace(0, log10(500), NL);

orderlist = {
    {[3, 0], [3, 0]};
    {[3, 1], [3, 1]}
    };
Norder = length(orderlist);


%analysis of proximal heavy ball algorithm
parfor i = 1:NL
    m= 1;
    L = Llist(i);

    
    %different operator classes for \nabla f
    op1_sml = op_sml(m, L);
    op1_quad = op_quad(m, L);
    
    %indicator function for \partial g
    op2 = op_sml(0, inf, 1);

    %proximal heavy ball
    gamma = 2/(m+L);
    lambda = 0.65;
    
    A = [1+lambda, -lambda; 1, 0];
    B = [-gamma, -gamma; 0 , 0];
    C = [1, 0; 1+lambda, -lambda];
    D = [0, 0; -gamma, -gamma];
    
    sK = ss(A, B, C, D, 1);

    %form the system
    sys_sml = opt_system({op1_sml,  op2}, [], sK);
    sys_quad = opt_system({op1_quad, op2}, [], sK);

    %pose and solve the problem
    man_sml = opt_analysis(sys_sml);
    man_quad = opt_analysis(sys_quad);

    for j = 1:Norder   
        sol_sml = man_sml.bisect(orderlist{j});
        rho_sml(i, j) = sol_sml.rho;
    
        sol_quad = man_quad.bisect(orderlist{j});
        rho_quad(i, j) = sol_quad.rho;    
    end

end

%% plot the result
figure(4)
clf
hold on
plot(Llist, rho_sml, 'linewidth', 2)
plot(Llist, rho_quad, 'linewidth', 2)

plot([Llist(1), Llist(end)], [1, 1], ':', 'linewidth', 2, 'color', 0.5*[1,1,1])
xlim([Llist(1), Llist(end)]);

xlabel('$L$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\rho$', 'interpreter', 'latex', 'fontsize', 16)
set(gca, 'xscale', 'log')
lname = {'Order [3, 0], $S_{m, L}$', 'Order [3, 1], $S_{m, L}$', ...
    'Order [3, 0], Quadratic', 'Order [3, 1], Quadratic', '$\rho$=1'};
legend(lname, 'location','southeast', 'interpreter', ...
    'latex', 'fontsize', 16)