function [fig] = plot_algsim(s, obj_rec, zzrec)
%PLOT_ALGSIM Summary of this function goes here
%   Detailed explanation goes here
fig = figure(1);
% fig = figure;
clf

x0 = s.t(1);
xe = s.t(end);
% xlim([x0, xe

        FS = 16;
        tiledlayout(2, 3)
nexttile
plot(s.t, s.z')
xlabel('$k$', 'Interpreter','latex')
ylabel('$z$', 'Interpreter','latex')
title('Iterate', 'FontSize', FS, 'Interpreter','latex')
xlim([x0, xe]);

nexttile
if nargin==3
    semilogy(s.t, sqrt(sum((s.z_avg - zzrec).^2, 1)))
    xlabel('$k$', 'Interpreter','latex')
    ylabel('$|\bar{z} - z^*|$', 'Interpreter','latex')
    title('Distance to Optima', 'FontSize', FS, 'Interpreter','latex')
else
    plot(s.t, s.w')
    xlabel('$k$', 'Interpreter','latex')
    ylabel('$w$', 'Interpreter','latex')
    title('Subdifferential', 'FontSize', FS, 'Interpreter','latex')
end
xlim([x0, xe]);

nexttile
plot(s.t, s.x(:, 1:end-1)')
xlabel('$k$', 'Interpreter','latex')
ylabel('$x$', 'Interpreter','latex')
title('State', 'FontSize', FS, 'Interpreter','latex')
xlim([x0, xe]);

% semilogy(s.t, sqrt(sum((s.x - xstar).^2, 1))')
% xlabel('$k$', 'Interpreter','latex')
% ylabel('$|x - x^*|$', 'Interpreter','latex')
% title('State Distance', 'FontSize', FS, 'Interpreter','latex')

nexttile
semilogy(s.t, abs(s.z_consensus))
xlabel('$k$', 'Interpreter','latex')
ylabel('$|z - z_{avg}|$', 'Interpreter','latex')
title('Consensus Error', 'FontSize', FS, 'Interpreter','latex')
xlim([x0, xe]);

nexttile
semilogy(s.t, abs(s.optimality))
xlabel('$k$', 'Interpreter','latex')
ylabel('$|\mathbf{1}^\top w|$', 'Interpreter','latex')
title('Optimality', 'FontSize', FS, 'Interpreter','latex')
xlim([x0, xe]);

nexttile

if isfield(s, 'payoff')
    plot(s.t, s.payoff)
    title('Payoff', 'FontSize', FS, 'Interpreter','latex')
    ylabel('$U$', 'Interpreter','latex')
else
    if nargin < 2
        plot(s.t, s.f)
        ylabel('$f$', 'Interpreter','latex')
    else
        semilogy(s.t, abs(s.f - obj_rec))
        ylabel('$f - f^*$', 'Interpreter','latex')
    end
    title('Function Value', 'FontSize', FS, 'Interpreter','latex')
    
end
xlabel('$k$', 'Interpreter','latex')
xlim([x0, xe]);

end

