function [regulator] = form_internal_model(p)
%FORM_INTERNAL_MODEL create the internal model by solving the regulator
%equation. 
% 
%
%
%TODO: implement for repeated nonlinearities

osign = p.opts.onesign;



%% break down the plant structure
[A, B, C, D] = ssdata(p.P);
B1 = B(:, p.iw);
B2 = B(:, p.iu);
C1 = C(p.iz, :);
C2 = C(p.iy, :);
D11 = D(p.iz, p.iw);
D12 = D(p.iz, p.iu);
D21 = D(p.iy, p.iw);
D22 = D(p.iy, p.iu);

n = size(A, 1);

s = length(p.iz);
N = [eye(s-1); -ones(1, s-1)];

if isfield(p.opts, 'tracking') && ~isempty(p.opts.tracking)
    S = blkdiag(p.opts.tracking.Sbeta, eye(s-1));
    Rbeta = ones(s, 1)*p.opts.tracking.Rbeta;
    TRACK = 1;
else
    S = eye(s);
    Rbeta = ones(s, 1);
    TRACK = 0;
end

ns = size(S,1);
nw = ns - s + 1;

%% solve the regulator equation for constant dynamics (static optimization)



%TODO: constant dynamics
if TRACK
    %dynamical tracking, complicated Rbetalvester
    reg_ans_mat = [zeros(n, nw), -B1*N; osign * Rbeta -D11*N];
    
    reg_ans = reshape(reg_ans_mat, [], 1);
    
    reg_mat_L = [A, B2; C1, D12];
    reg_mat_RR = S;
    reg_mat_RL = blkdiag(speye(size(A, 1)), sparse(s, s));


    reg_mat = kron(speye(ns), reg_mat_L) - kron(reg_mat_RR', reg_mat_RL);
   
    try 
        reg_sol_vec = reg_mat \ reg_ans;
        % reg_sol = 
        % reg_Rbetalv = Rbetalvester(Rbetalv_A, Rbetalv_B, Rbetalv_C);
        % Rbetalv_A = 
        % reg_sol = Rbetalvester()      
        reg_err = norm(reg_mat * reg_sol_vec - reg_ans);
        reg_sol = reshape(reg_sol_vec, size(reg_mat_L, 1), []);
    catch
        warning('Regulator equation cannot be solved')
    end

    if isnan(reg_err) || norm(reg_err) > 1e-7
        warning('Regulator equation cannot be solved')
        return
    end

    Pi = reg_sol(1:n, :);
    Gam = reg_sol(n+1:end, :);
    
    Phi = D21 * [zeros(s, nw), N] + D22*Gam + C2*Pi;

    reg_err_mat = reg_mat_RL * [Pi; Gam]*S - reg_mat_L * [Pi; Gam] + reg_ans_mat;
else
    %constant tracking, simple Rbetalvester
    if n == 0
        Pi = [];
        % Gam = D12 \ [-ones(s, 1), -D11*N];
        Gam = D12 \ [osign * ones(s, 1), -D11*N];
        Phi = D21 * [zeros(s, 1), N] + D22*Gam;
    else
        %TODO: THE ONES SIGN
        reg_ans = [zeros(n, 1), -B1*N; osign * ones(s, 1), -D11*N];
        reg_mat = [A - eye(n), B2; C1, D12];
        
        try 
            reg_sol = reg_mat \ reg_ans;
            % TODO: replace by Rbetalvester. This is an explicit solution for
            % static optimization
    
            reg_err = norm(reg_mat * reg_sol - reg_ans);
        catch
            warning('Regulator equation cannot be solved')
        end
    
        if isnan(reg_err) || norm(reg_err) > 1e-7
            warning('Regulator equation cannot be solved')
            return
        end
    
        %recover basic solution to the regulator equation
        Pi0 = reg_sol(1:n, :);
        Gam0 = reg_sol(n+1:end, :);        
        Phi0 = D21 * [zeros(s, 1), N] + D22*Gam0 + C2*Pi0;


        if isfield(p.opts, 'basis_trans') && p.opts.basis_trans
            
            %Apply coordinate transformation
            %to get rid of feedback term in P
            P1 = null(Phi0);
            P2 = orth(Phi0);
            P = [P1, P2];
    
            Pi = Pi0 * P;
            Gam = Gam0 * P;
            Phi = Phi0 * P;
            
        else
            Pi = Pi0;
            Gam = Gam0;
            Phi = Phi0;
            
        end

    end
end
%% build the internal model


Abar = S;
Bblock = [eye(ns), zeros(ns, s), zeros(ns, s)];
Cblock = [-Gam;
          Phi];

Dblock = [zeros(s, ns), eye(s), zeros(s);
    zeros(s, ns), zeros(s), eye(s)];


model = ss(Abar, Bblock, Cblock, Dblock, 1);

%label all attributes
for i = 1:nw - 1+1
    model.StateName{i} = sprintf('shift%d', i);
end
for i = nw - 1+2:ns
    model.StateName{i} = sprintf('subgrad%d', i);    
end  
 

for i = 1:(3*s + nw - 1)    
    if i<=s + nw - 1   
        model.InputName{i} = sprintf('u1_%d', i); %sprintf('u1_%d', i);
    elseif i <= 2*s + nw - 1
        model.InputName{i} = sprintf('u2_%d', i - s);    
    else            
        model.InputName{i} = sprintf('w%d', i- 2*s);
    end
end

for i = 1:(2*s)
    if i <= s
        model.OutputName{i} = sprintf('z%d', i);
    else
        model.OutputName{i} = sprintf('y%d', i-(s));
    end
end

%internal model
model = model(:, [(2*s+nw - 1+1):(3*s+nw - 1), 1:(2*s+nw - 1)]);

%interconnected with Rbetastem
model_sys = lft(p.P, model, s, s);

for i = 1:n
    model_sys.StateName{i} = sprintf('x%d', i);
end



%% package and export
regulator = struct('S',S, 'Pi', Pi, 'Gam', Gam, ...
    'Phi', Phi, 'model', model, 'model_sys', model_sys);

end

