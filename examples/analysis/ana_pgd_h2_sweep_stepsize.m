%define the operators
m = 1; L = 10;
op1= op_pcc(); %pcc goes first, to not take information from the noisy gradient
op2= op_sml(m, L);
ops = {op1, op2};


%gain from (error in w2) to (z2 - z*)
network = bridge_pass_through(2);
[network, iwp] = network.add_oracle_input(2, []);
[network, izp] = network.perf_output_z(2);

%define specification
GAIN = 10;
spec_stoch = spec_h2(GAIN, 1, iwp, izp);
spec_stoch.target = true;
specs = {spec_stoch};


%sweep over the pgd stepsize
gam_max = 2/(L + m);
order = [1,1];

Ngam = 400;
gsweep = linspace(0, gam_max, Ngam+1);
gsweep = gsweep(2:end);

h2gain = zeros(Ngam, 1);
rho = zeros(Ngam, 1);

h2gain_joint = zeros(Ngam, 1);
rho_joint = zeros(Ngam, 1);

parfor i = 1:Ngam
    %PGD algorithm
    gamma = gsweep(i);
    K = ss(1, [-gamma, -gamma], [1; 1], [-gamma, 0; -gamma, 0],1);      
    
    %form the systems
    sys_clean = opt_system(ops, [], K);
    sys_noisy = opt_system(ops, network, K);
    
    %pose the analysis problems

    %for convergence rate
    man_noisy= opt_analysis(sys_noisy);
    sol_noisy = man_noisy.solve_single(order, specs);
    
    %for the sensitivity
    man_clean = opt_analysis(sys_clean);
    sol_clean = man_clean.bisect(order);


    %bisection on convergence rate, inner target sensitivity
    specs_joint = {spec_stability(1), spec_stoch};
    man_joint= opt_analysis(sys_noisy);
    sol_joint= man_noisy.bisect(order, specs_joint);
    
    
    h2gain(i) = sqrt(sol_noisy.objective);
    rho(i) = sol_clean.rho;

    h2gain_joint(i) = sqrt(sol_joint.objective);
    rho_joint(i) = sol_joint.rho;
end

save("ana_pgd_h2_sweep")

%% plot the result (pareto)

figure(3)
clf
hold on

plot(h2gain, rho, 'LineWidth',2)
plot(h2gain_joint, rho_joint, 'LineWidth',2)

xlabel('Sensitivity', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\rho$', 'interpreter', 'latex', 'fontsize', 16)
rho_best = (L-m)/(L+m);
xl = xlim;
plot(xl, rho_best * [1, 1], ':', 'color', 0.5 * [1, 1, 1], 'LineWidth', 2);
lname = {'Separate Search', 'Joint Search', '$\rho = \frac{L-m}{L+m}$'};
legend(lname, 'location','northeast', 'interpreter', ...
    'latex', 'fontsize', 16)

text(0.15, rho_best+0.02, '$\gamma = \frac{2}{m+L}$', 'Interpreter','latex', 'fontsize', 14)
text(0+0.03, 0.98, '$\gamma = 0$','Interpreter','latex', 'fontsize', 14)

%% plot v.s. gamma
figure(4)
clf
% tiledlayout(2, 1)
% nexttile
hold on
plot(gsweep, rho, 'linewidth', 2)

xlim([0, gam_max])

xlabel('Stepsize $\gamma$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\rho$', 'interpreter', 'latex', 'fontsize', 16)


yyaxis right
plot(gsweep, h2gain, '--', 'linewidth', 2)
ylabel('$\ell_2$ gain', 'interpreter', 'latex', 'fontsize', 16)

set(gca, 'xscale', 'log')
lname = {'$\rho$', '$\ell_2$ gain'};
legend(lname, 'location','southeast', 'interpreter', ...
    'latex', 'fontsize', 16)

% nexttile 
% plot(gsweep, rho)


%create parto