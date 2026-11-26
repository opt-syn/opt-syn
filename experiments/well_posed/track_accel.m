%track an accelerating object
%identify the internal model structure


%% exosystem (accelerating optimal point)
%order of the system
% r = 4;
r = 2;
dt = 1;
% dt = 0.5;

Sbeta = zeros(r);
for i = 1:r
    for j = 1:r
        if j >= i
            Sbeta(i, j) = dt^(j-i);
        end
        if j > i
            Sbeta(i, j) = Sbeta(i, j) / (j-i);
        end
    end
end

Rbeta = zeros(1, r);
Rbeta(1) = 1;

%% algorithm
m = 1;
L = 10;
op = op_sml(m, L);

gamma = 2/(m+L);
K = ss(1, -gamma, 1, 0);

sys = opt_system(op, [], K);
sys.tracking = struct('Sbeta', Sbeta, 'Rbeta', Rbeta);



%% regulator check
reg = regulator_lti(sys);
regcl = reg.check_regulator()