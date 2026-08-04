%pose the operator class
m = 1; L = 5;
op1 = op_sml(m, L);

op2 = op_pcc();
ops = {op1, op2};


%use a transfer function representation 
%to model the channel memory
z = tf('z', 1);

%sweep over alpha and order
Nalpha = 100;
alpha_list = logspace(-2, 1, Nalpha);
orderlist = {0, [1, 0], [0, 1], [1, 1], [2, 0], [0,2]};
Norder = length(orderlist);
rho_list = zeros(Nalpha, Norder);


parfor i = 1:Nalpha
    alpha = alpha_list(i);
    ascale = (2*(alpha+1));
    
    %form the channel memory network
    P = [0, 0, z/(z+alpha), 0;
        0 , 0, 0, 1;
        z/(z+alpha), 0, 0, 0;
        0, 1, 0, 0];
    
    network = genplant(P);
    network.nw = 2; network.nu = 2;
    network.nz = 2; network.ny = 2;

    %create the controller
    gamma = 0.4;
    lambda = 0.2;

    K = ss([1], [-gamma*lambda, -gamma*lambda/(alpha+1)], ...
        [1+alpha; 1], [0, 0; -gamma, -gamma/(alpha+1)],1);
    
    sys = opt_system(ops, network, K);
    
    %solve Analysis at each order
    man = opt_analysis(sys);    
    for j = 1:Norder
        order = {orderlist{j}, orderlist{j}};
        sol = man.bisect(order);
        rho_list(i, j) = sol.rho;
    end
       
end


%% plot the result
figure(5)
clf
cc = linspecer(Norder+1);
hold on


% j = Norder;
jlist = [1, Norder];
for i = 1:length(jlist)
    plot(alpha_list, rho_list(:, jlist(i)), 'linewidth', 2, 'color', cc(1+i, :))
end
plot(alpha_list, ones(size(alpha_list)), ':', 'LineWidth', 2, 'color', cc(1, :))
set(gca, 'Xscale', 'Log')
ylabel('$\rho$', 'interpreter', 'latex','Fontsize', 14)
xlabel('$\alpha$', 'interpreter', 'latex','Fontsize', 14)

legend({'Order = [0,0]', 'Order = [0,2]', '$\rho=1$'}, 'location','northwest', 'interpreter', 'latex','Fontsize', 14)
% legend({'Order [1, 1]', '$\rho=1$'}, 'location','northeast', 'interpreter', 'latex')
