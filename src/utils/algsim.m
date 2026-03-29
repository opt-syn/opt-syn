function sim=algsim(p);
%function sim=algsim(p);
%
%Simulate closed loop system of an algorithm.
%
%This assumes Kronecker structure A = A_a kron I_d
%   TODO: handle other coordinate partitioning structures
%
%The theory is exposed in C.W. Scherer, Ch. Ebenbauer, T. Holicki, 
%Optimization Algorithm Synthesis based on Integral Quadratic Constraints: A Tutorial, 
%62nd IEEE Conference on Decision and Control. 
%
%An extended version is available on arXiv under https://doi.org/10.48550/arXiv.2306.00565
%All references in the code are related to these paper.
%
%Generate trajectories through recursion (Figure 1)
%
%x_(t+1)=Ax_t+Bw_t, w_t^i=partial f_i(z_t^i)
%z_t    =Cx_t+Dw_t 
%
%with initial condition zero. A nonzero trajector is 
%generated with a function f for which nabla f(0)<>0.
%
%The input is collected in the structure p with following fields:
%p.alg            Closed loop system in ss format as used for simulation 
%p.d              Dimension for Kroneckering 
%p.grad           Matlab function to compute gradient 
%
%%Results are collected in structure s with following additional fields:
%sim.x              state trajectory
%sim.w              input trajectory
%sim.z              output trajectory
%sim.t              time instants

[A,B,C,D]=ssdata(p.alg);
% if norm(d(1,1))>1e-6;
%     error('Direct feedthrough not zero.');
% end;

% A=kron(a,eye(p.d));
% B=kron(b,eye(p.d));
% C=kron(c,eye(p.d));
% D=kron(d,eye(p.d));
% ND = size(d, 1);

s = size(D, 1);


n=size(A,1);

if isfield(p,'x');
    if length(p.x) == n*p.d
        x_vec=p.x;
    else
        x_vec=kron(p.x,ones(p.d,1));
    end
else
    x_vec=zeros(n*p.d,1);
end

x = reshape(x_vec, [p.d, n])';

wp = zeros(s, p.d);
zp = zeros(s, p.d);
fp = zeros(s, p.T);
z_avg = zeros(p.d, p.T);
z_con = zeros(p.d*s, p.T);
w_opt = zeros(p.d, p.T);
w_opt_err = zeros(1, p.T);
z_con_err = zeros(1, p.T);
if isfield(p, 'payoff')
    pfv = [];
    PAYOFF = 1;
else
    PAYOFF = 0;
    pfv = [];
end

xv = [x(:)];
wv = [];
zv = [];
fv = [];
for t=1:p.T
    for i = 1:s
        %ASSERT D is lower-triangular
        if i==1
            %flush the w and z values
            wp = NaN*wp;
            zp = NaN*zp;
            fp = NaN*fp;
        end
        % dind_x = (1:p.d) + (i-1)*p.d;
        
        
        vi = C(i, :)*x;
        if i > 1
            % dind_w = (1:(i-1)*p.d); %there was a bug here, 
            % the entire D matrix up to i-1 was not used.
            % vwi = D(dind_x, dind_w)*wp(dind_w);
            vwi = D(i, 1:(i-1)) * wp(1:(i-1), :);
            vi = vi + vwi;
        end

        %prox or grad?
        % vi_vec = reshape(vi', [], 1);
        if D(i, i) ~= 0
            

            zi = p.prox{i}(t, -D(i, i), vi')';   
            
            wi = (vi - zi)/(-D(i, i));
            % wi = reshape(wi_vec, [], s)';
        else
            zi = vi;
            % zi_vec = reshape(vi', [], 1);
            wi= p.grad{i}(t, vi');
            % wi = reshape(wi_vec, [], s)';
        end

        if PAYOFF 
            pfp = p.payoff(t, zi');
        else
            pfp = [];
            fp(i, t) = p.f{i}(t, zi');
        end
        
        wp(i, :) = wi;
        zp(i, :) = zi;
    end
    % zv=[C*x+D*w];
    % w=p.grad(C(1:p.d,:)*xv(:,end));
    xp=A*x+B*wp;
    % zp=C*xv(:,end)+D*w;
    % wp=p.grad(zp(1:p.d));
    
    xv=[xv xp(:)];
    zv=[zv zp(:)];
    wv=[wv wp(:)];   
    
    
    z_avg(:, t) = sum(zp, 1)' * (1/s);
    z_con_curr = zp - ones(s, 1) * z_avg(:, t)';
    z_con(:, t) = reshape(z_con_curr', [], 1);
    w_opt(:, t) = sum(wp, 1)';

    z_con_err(t) = norm(zp - ones(s, 1) * z_avg(:, t)' , 'fro');
    w_opt_err(t) = norm(w_opt(:, t), 'fro');
    fv =[fv sum(fp(:, t))];
    pfv = [pfv pfp];
    x = xp;
end;
sol = struct;

sim.x=xv;
sim.w=wv;
sim.z=zv;
sim.t=0:p.T-1;
sim.f = fv;
if PAYOFF
    sim.payoff = pfv;
end
%consensus
sim.z_avg = z_avg;
sim.z_consensus = z_con;
sim.optimality = w_opt;
sim.z_consensus_error = z_con_err;
sim.optimality_error = w_opt_err;