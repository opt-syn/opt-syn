%define the operators
m = [1, 1, 1, 1];
L = [2, 4, 6, 8];
s = length(m);

op_list = cell(s, 1);
for i = 1:s
    op_list{i}= op_sml(m(i), L(i));       
end

%parameters of the controller, independent of omega
b0 = [-0.05; -0.1; -0.05];
b1 = -0.05;
b2 = -0.1;

BK = b0 * ones(1, s);
CK = ones(s, 1) * Rbeta;

DK = zeros(s);
DK(end, :)= DK(end, :) + b1;
DK(:, 1)= DK(:, 1) + b2;


%sweep over the angle
Nomega = 320;
omegalist = pi * linspace(-1, 1, Nomega);

rholist = zeros(1, Nomega);

order = {2, 2, 2, 2};
parfor i = 1:length(omegalist)

    %define the path of the optimal solution
    theta = omegalist(i);
    Sbeta = blkdiag(1, givens(cos(theta), sin(theta)));
    Rbeta = [1, 1, 0];

    
    %create the controller    
    K = ss(Sbeta, BK, CK, DK, 1);
    
    
    %form the system
    sys = opt_system(op_list, [], K);
    
    
    %run Analysis and store result
    man = opt_analysis(sys);  
    sol = man.bisect(order);

    rholist(i) = sol.rho;    
end

%% plot the result


%save prior runs, e.g. rholist_0_0
figure(1)
clf
hold on
cc = linspecer(4);
plot(omegalist, rholist_0_0, 'linewidth', 2, 'color', cc(1, :))
plot(omegalist, rholist_1_0, 'linewidth', 2, 'color', cc(2, :))
plot(omegalist, rholist_2_0, 'linewidth', 2, 'color', cc(3, :))
plot([-pi, pi], [1, 1], ':', 'linewidth', 2, 'color', 0.5*[1,1,1])
xlabel('$\omega$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\rho$', 'interpreter', 'latex', 'fontsize', 16)
xlim([min(omegalist), max(omegalist)]);
legend({'Order [0,0]', 'Order [1, 0]', 'Order [2, 0]', '$\rho$=1'}, 'location','northeast', 'interpreter', ...
    'latex', 'fontsize', 16,'location', 'north')


%% compare against quad
figure(2)

plot(omegalist, rholist_0_0, 'linewidth', 2, 'color', cc(1, :))
plot(omegalist, rholist-rholist_0_0, 'linewidth', 2, 'color', cc(2, :))
