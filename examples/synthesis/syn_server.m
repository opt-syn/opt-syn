%form the operators
m = [0, 1, -2, 1, 1, 0];
L = [5, 2, 1, 6, 1, inf];
s = length(m);

ops= cell(s, 1);
for i = 1:s
    ops{i} = op_sml(m(i), L(i));
end

%form the systems
sys = opt_system(ops);

%solve the problem
man = opt_synthesis(sys);
sol = man.bisect();  %0.7679

%% plot solutions