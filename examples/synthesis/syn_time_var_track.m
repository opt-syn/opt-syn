%tracking a rotating object with time-varying delay
m = [1];
L = [1.5];
s = length(m);

TRACK = true;

if TRACK

    theta = pi/16;
    Sbeta = blkdiag(1, givens(cos(theta), sin(theta)));
    Rbeta = [1, 1, 0];
else
    Sbeta = 1;
    Rbeta = 1;

end
op_list = cell(s, 1);

for i = 1:s
    op_list{i}= op_sml(m(i), L(i));   
end
 

%delay network
delay_max = 1;
[Pprim, Gcon] = delay_primitives(0:delay_max, 0, -1:0:1, 0);
% [Pprim, Gcon] = delay_primitives(0:delay_max, 0:delay_max, -1:0:1, -1:0:1);
% 
% Pdelay = cell(s, 1);
% for i = 1:delay_max+1
%     for j = 1:s
%         if j==1
%             Pdelay{i} = Pprim{i};
%         else
%             Pdelay{i} = blkdiag(Pdelay{i}, Pprim{i});
%         end
%     end
% end
network = genplant_poly(Pprim);
% network = genplant_poly(Pdelay);
Gsnap = delay_snap_graph(delay_max, 0);
% Gsnap  = ones(delay_max + 1);
% Gsnap = eye(delay_max);

% Gsnap = eye(delay_max+1);
%form the system
sys = opt_system_switched(op_list, network, [], Gsnap);
sys.tracking = struct('Sbeta', Sbeta, 'Rbeta', Rbeta);

reg = regulator_switched(sys);

man = opt_synthesis(sys);

spec = spec_stability(0.99);
sol = man.solve_single([], {spec});
sol.rho


%% plot the result
if ~sol.status
d = 10;
Q = rand_quad(d, m, L);

%tracking the solution
if TRACK
    SA = kron(Sbeta, eye(d));
    SAy = kron(Rbeta, eye(d));
    eta0 = [randi(201, [d, 1]) - 100; 
        (randi(31, [2*d, 1]) + 60).*sign(2*rand(2*d, 1) - 1)];
    
    shift = @(k) SAy * (SA^k) * eta0;

else
    shift = @(k) 0;
end
bstar_center = randi(101, [d, 1]) - 50;
bstar = @(k) bstar_center + shift(k);
op1 = op_sim_quad(Q, bstar);


sys_snap = sol.export_sim({op1});
% sys_snap.adj = eye(2);
sim_snap = alg_sim(sys_snap, d);
T = 200;
sim_out_snap = sim_snap.sim(T);


% plot
plt_snap = alg_plotter(sim_out_snap);
plt_snap.plot({'xn', 'w', 'res_w', 'xc', 'z', 'delay'}, 12)
end