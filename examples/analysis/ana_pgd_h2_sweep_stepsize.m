%define the operators
m = 1;L = 10;
op1= op_pcc(); %pcc goes first, to not take information from the noisy gradient
op2= op_sml(m, L);
ops = {op1, op2};


%gain from (error in w2) to (z2 - z*)
network = bridge_pass_through(2);
[network, iwp] = network.add_oracle_input(2, []);
[network, izp] = network.perf_output_z(2);

%define specification
GAIN = 10;
spec = spec_h2(GAIN, 1, iwp, izp);
spec.target = true;
specs = {spec};


%sweep over the pgd stepsize
gam_max = 2/(L + m);
order = [1,1];


Ngam = 200;
gsweep = linspace(0, gam_max, Ngam+1);
gsweep = gsweep(2:end);

l2gain = zeros(Ngam, 1);
rho = zeros(Ngam, 1);

parfor i = 1:Ngam

    gamma = gsweep(i);
    %PGD algorithm
    K = ss(1, [-gamma, -gamma], [1; 1], [-gamma, 0; -gamma, 0],1);      
    
    %form the systems
    sys_clean = opt_system(ops, [], K);
    sys_noisy = opt_system(ops, network, K);
    
    %pose the analysis problem

    
    man_noisy= opt_analysis(sys_noisy);
    sol_noisy = man_noisy.solve_single(order, specs);
    
    man_clean = opt_analysis(sys_clean);
    sol_clean = man_clean.bisect(order);
    
    
    l2gain(i) = sqrt(sol_noisy.objective);
    rho(i) = sol_clean.rho;
end

%% plot the result

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
plot(gsweep, l2gain, '--', 'linewidth', 2)
ylabel('$\ell_2$ gain', 'interpreter', 'latex', 'fontsize', 16)

set(gca, 'xscale', 'log')
lname = {'$\rho$', '$\ell_2$ gain'};
legend(lname, 'location','southeast', 'interpreter', ...
    'latex', 'fontsize', 16)

% nexttile 
% plot(gsweep, rho)