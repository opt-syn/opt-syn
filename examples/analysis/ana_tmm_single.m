%define the operators
m = 1;L = 10;
op1= op_sml(m, L);
op2 = op_sml(0, inf);
ops = {op1, op2};

%define the controller (unconstrained and constrained/composite)
K_uncon = tmm(m, L);
K_con = tmm_prox(m, L);

%form the system
sys_uncon = opt_system(op1, [], K_uncon);
sys_con= opt_system(ops, [], K_con);

%pose the analysis problem
order = [1,1];

man_uncon = opt_analysis(sys_uncon);
man_con = opt_analysis(sys_con);

%solve the analysis problem
sol_uncon = man_uncon.bisect(order);
sol_con = man_con.bisect(order);

best_rho = 1 - sqrt(m/L);
[best_rho, sol_uncon.rho, sol_con.rho]
%0.6838    0.6838    0.6838
 

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
    %composite triple momentum
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