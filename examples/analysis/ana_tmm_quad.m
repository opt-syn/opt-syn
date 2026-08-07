
%proximal triple-momentum method from https://arxiv.org/abs/2605.22929
%both have the convergence rate (1 - sqrt(m/L))

%define the operators
m = 1;
L = 10;

%different operator classes for op1
op1= op_quad(m, L);
op2 = op_sml(0, inf);
ops = {op1, op2};

%define the controller (composite)


%form the system
K_uncon = tmm(m, L);
K_con = tmm_prox(m, L);
sys_uncon = opt_system(op1, [], K_uncon);
sys_con= opt_system(ops, [], K_con);


order = [1, 1];

man_uncon = opt_analysis(sys_uncon);
man_con = opt_analysis(sys_con);

sol_uncon = man_uncon.bisect(order);
sol_con = man_con.bisect(order);

best_rho = 1 - sqrt(m/L);
[best_rho, sol_uncon.rho, sol_con.rho]
%order vs. rho:
%[0, 0]: 0.6838    1.0682    1.1022
%[1, 0]: 0.6838    0.6838    0.7277
%[0, 1]: 0.6838    1.0682    1.0682
%[1, 1]: 0.6838    0.6838    0.6840
%[2, 0]: 0.6838    0.6838    0.7277
%[2, 1]: 0.6838    0.6838    0.6838
%order >= [1, 1] seems to give best results


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