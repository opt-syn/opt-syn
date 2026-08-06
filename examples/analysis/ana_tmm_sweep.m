%composite triple-momentum method from https://arxiv.org/abs/2605.22929
%both have the convergence rate (1 - sqrt(m/L))


Nm = 300;
mlist = linspace(0, 1, Nm+1);
mlist = mlist(2:end);
rho = zeros(Nm, 2, 3);

parfor i = 1:Nm
    for j = 1:2
        if j==1
            m = 1;
            L = 1/mlist(i);
        else
            m = mlist(i);
            L = 1;
        end
        
        %define the operators
        op1= op_sml(m, L);
        op2 = op_sml(0, inf);
        ops = {op1, op2};
        
        %define the controller (composite)              
        K_uncon = tmm(m, L);
        K_con = tmm_prox(m, L);
        
        %form the system
        sys_uncon = opt_system(op1, [], K_uncon);
        sys_con= opt_system(ops, [], K_con);
        
        
        %pose and solve the problem
        order = [1,1];
        
        man_uncon = opt_analysis(sys_uncon);
        man_con = opt_analysis(sys_con);
        
        sol_uncon = man_uncon.bisect(order);
        sol_con = man_con.bisect(order);
    
        best_rho = 1 - sqrt(m/L);
    
        rho(i,j,  :) = [sol_uncon.rho, sol_con.rho, best_rho];
    end
end

%% plot the sweep
figure(1)
clf
tiledlayout(2, 2)
for j = 1:2
nexttile
hold on
if j == 1
    xla = "$m$";
else
    xla = "$1/L$";
end
plot(mlist, squeeze(rho(:, j,  1:2)), 'linewidth', 2)
plot(mlist, squeeze(rho(:, j, 3)), '--', 'linewidth', 2, 'color', 0.5*[1,1,1]);
xlabel(xla, 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\rho$', 'interpreter', 'latex', 'fontsize', 16)
xlim([0, 1]);
ylim([0, 1])
if j ==1
legend({'TMM', 'Composite TMM', '$1-\sqrt{m/L}$'}, 'location','northeast', 'interpreter', ...
    'latex', 'fontsize', 16,'location', 'northeast')
end
nexttile
hold on
plot(mlist, squeeze(rho(:, j, 1:2)-rho(:, j, 3)), 'linewidth', 2)
xlabel(xla, 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\rho - \left(1-\sqrt{m/L}\right)$', 'interpreter', 'latex', 'fontsize', 16)
xlim([0, 1]);
set(gca, 'yscale', 'log')
end

%% functions to generate algorithms
function K = tmm(m, L)
    %triple momentum algorithm
    rho=1-1/sqrt(L/m);
    al=(1+rho)/L;
    be=rho^2/(2-rho);
    ga=rho^2/((1+rho)*(2-rho));
    A=[1+be -be;1 0];B=[-al;0];C=[1+ga -ga];
    K = ss(A,B,C,0,1);
end

function K = tmm_prox(m, L)
    %triple momentum proximal algorithm
    A = [(1.0 - sqrt(m / L)) / (1.0 + sqrt(m / L)),  2.0 * sqrt(m / L) / (1.0 + sqrt(m / L));
        sqrt(m / L) * (1.0 - sqrt(m / L)) / (1.0 + sqrt(m / L)), 1.0 - sqrt(m / L) + 2.0 * m / L / (1.0 + sqrt(m / L))];

    B = [-1.0 / L, -1.0 / L;
        -1.0 / (sqrt(m / L) * L), -1.0 / (sqrt(m / L) * L)];
    C = [    (1.0 - sqrt(m / L)) / (1.0 + sqrt(m / L)),2.0 * sqrt(m / L) / (1.0 + sqrt(m / L));   
        sqrt(m / L) * (1.0 - sqrt(m / L)) / (1.0 + sqrt(m / L)),    1.0 - sqrt(m / L) + 2.0 * m / L / (1.0 + sqrt(m / L))];

    D = [0.0, 0.0;
        -1.0 / (sqrt(m / L) * L), -1.0 / (sqrt(m / L) * L)];

    K = ss(A, B, C, D, 1);
end