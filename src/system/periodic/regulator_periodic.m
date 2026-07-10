classdef regulator_periodic < regulator_switched
    %REGULATOR_PERIODIC Regulator for periodic systems    
    %
    % [x(k+1)] = [A(k)    Bd(k)    Bu(k)  ][x(k)]   state transition
    % [e(k)  ] = [Ce(k)   Ded(k)   Deu(k) ][d(k)]   output to  regulated error    
    % [zp(k) ] = [Cy(k)   Dyd(k)   Dyu(k) ][u(k)]   output to controller    

    %A(k) = A(k+T) for some known time T
    %
    %instances of these algorithms include cyclic coordinate descent
    %methods. Periodic systems can also be unrolled into an LTI system
    %(monodromy methods): a single large LMI system rather than multiple 
    % coupled smaller LMI systems

    methods
        function obj = regulator_periodic(sys)
            %REGULATOR_PERIODIC undefined
            %   undefined
            obj@regulator_switched(sys)
        end

    %     function ns = Nss(obj)
    %         %NSS: number of subsystems
    %         ns = obj.sys.Nss;
    %     end
    % 
    %     function obj = form_internal_model(obj)
    %         %FORM_INTERNAL_MODEL create the internal model by solving the regulator
    %         %equation. Inputs are the system (P, bind, tracking, op)
    %         %
    %         %op is important for which oracles are equaltiy constarints and
    %         %which are inequality constraints
    % 
    %         Npre = obj.sys.get_consensus(obj.sys.op, obj.sys.bind);
    %         c = obj.sys.op{1}.c; %coordinate lifts: change this later?
    %         N = kron(Npre, eye(c));
    % 
    %         [sN0, dN0] = size(Npre);
    %         [sN, dN] = size(N);
    %         n = obj.sys.P.nx;
    %         nu = obj.sys.P.nu;
    % 
    %         [Sbeta, Rbeta] = obj.sys.get_tracked_opt();
    % 
    %         if isempty(obj.sys.tracking)
    %             S = eye(size(N, 2)+c);
    %             R = S;
    %             reg_ans = [];
    %             reg_mat = [];
    % 
    %             %go through each subsystem
    %             for i = 1:obj.Nss
    %                 Pcurr = obj.sys.P{i};
    %                 [A, B1, B2, C1, D11, D12, C2, D21, D22] = Pcurr.ss_zy_wu();
    % 
    %                 reg_ans_curr = [zeros(n, c), -B1*N;  -kron(ones(sN0,1), eye(c)), -D11*N];
    %                 reg_mat_curr = [A, B2; C1, D12];
    % 
    %                 reg_ans = [reg_ans; reg_ans_curr];
    %                 reg_mat = blkdiag(reg_mat, reg_mat_curr);
    %             end
    % 
    %             %now collect the shifts
    % 
    %             eye_shift =blkdiag(eye(n), zeros(size(D12)));
    %             next_trans = kron(eye(obj.Nss), eye_shift);
    % 
    %             next_trans = circshift(next_trans, size(eye_shift, 2), 2);
    % 
    %             reg_mat = reg_mat - next_trans;
    % 
    %             %solve the regulator equation
    %             null_basis = null(reg_mat, 'rational');
    %             % try                    
    %                 sol0 = lsqminnorm(reg_mat, reg_ans);
    %             % catch
    %             sol_err = reg_mat * sol0 - reg_ans;
    %             if norm(sol_err) > 1e-12
    %                 error('Regulator equation cannot be solved')
    %             end
    % 
    %             nnull = size(null_basis, 2);
    % 
    %             %extract the solution
    % 
    %             Pi0 = cell(obj.Nss, 1);
    %             Gam0 = cell(obj.Nss, 1);
    %             Phi0 = cell(obj.Nss, 1);
    % 
    %             Pi_basis = cell(obj.Nss, 1);
    %             Gam_basis = cell(obj.Nss, 1);
    %             Phi_basis = cell(obj.Nss, 1);
    % 
    %             count = 0;
    %             for i = 1:obj.Nss
    %                 %get the regulator equation solution
    %                 Pcurr = obj.sys.P{i};
    %                 [~, ~, ~, ~, ~, ~, C2, D21, D22] = Pcurr.ss_zy_wu();
    % 
    % 
    %                 ind_pi = count + (1:n);
    %                 ind_gam = count + n+ (1:nu);
    %                 Pi0{i} = sol0(ind_pi, :);
    %                 Gam0{i} = sol0(ind_gam, :);
    %                 Phi0{i} = D21 * [zeros(sN, c), N] + D22*Gam0{i} + C2*Pi0{i};
    % 
    %                 %get the free parameters
    %                 if nnull
    %                     Pi_basis_pre = null_basis(ind_pi, :);
    %                     Gam_basis_pre = null_basis(ind_gam, :);
    %                     Phi_basis_pre = D22*Gam_basis_pre + C2*Pi_basis_pre;
    % 
    %                     Pi_basis{i} = kron(Pi_basis_pre, eye(nnull));
    %                     Gam_basis{i} = kron(Gam_basis_pre, eye(nnull));
    %                     Phi_basis{i} = kron(Phi_basis_pre, eye(nnull));
    %                 else
    %                     Pi_basis{i} = [];
    %                     Gam_basis{i} = [];
    %                     Phi_basis{i} = [];
    %                 end
    % 
    %                 count = count + n + nu;
    % 
    %             end                                    
    % 
    % 
    %         else
    %             error('Periodic regulation: tracking not yet supported')
    %             %TODO: to be implemented
    %         end
    % 
    % 
    %         obj.S = S;
    %         obj.R = R;
    %         obj.Pi = Pi0;
    %         obj.Gam = Gam0;
    %         obj.Phi = Phi0;
    %         obj.Pi_basis = Pi_basis;
    %         obj.Gam_basis = Gam_basis;
    %         obj.Phi_basis = Phi_basis;
    %     end
    % 
    %     function [regulator_closed] = check_regulator(obj)
    %         %CHECK_REGULATOR is the regulator equation satisfied?
    %         sys_cl = lft(obj.sys.P, obj.sys.K);
    % 
    %         Npre = obj.sys.get_consensus(obj.sys.op, obj.sys.bind);
    %         c = size(sys_cl{1}.D, 1)/length(obj.sys.bind); %coordinate lifts: change this later?
    %         N = kron(Npre, eye(c));
    % 
    % 
    %         if isempty(obj.sys.tracking)
    %             S = eye(size(N, 2)+1);
    %             R = S;
    %             reg_ans = [];
    %             reg_mat = [];
    % 
    %             n = length(sys_cl{1}.A);
    % 
    %             %go through each subsystem
    %             for i = 1:obj.Nss
    %                 S = eye(1+size(N, 2));
    %                 R = S;
    % 
    %                 [A, B, C, D] = ssdata(sys_cl.P{i});
    % 
    % 
    %                 reg_ans = [zeros(n, 1), -B*N;  -ones(sN, 1), -D*N];
    %                 reg_mat = [A - eye(n); C];
    % 
    %                 reg_ans = blkdiag(reg_ans, reg_ans_curr);
    %                 reg_mat = blkdiag(reg_mat, reg_mat_curr);
    %             end
    % 
    %             %now collect the shifts                
    %             next_trans = kron(eye(obj.Nss), -eye(n));
    %             next_trans = circshift(next_trans, -n, 1);
    % 
    %             reg_mat = reg_mat - next_trans;
    % 
    %             %solve the regulator equation
    %             null_basis = null(reg_mat, 'rational');
    %             try                    
    %                 sol0 = reg_mat \ reg_ans;
    %             catch
    %                 warning('Regulator equation cannot be solved')
    %             end
    % 
    %             nnull = size(null_basis, 2);
    % 
    %             %extract the solution
    % 
    %             Pi0 = cell(obj.Nss, 1);
    %             Th0 = cell(obj.Nss, 1);
    % 
    % 
    % 
    %             count = 0;
    %             nxn = obj.sys.nxn;
    %             nxi = obj.sys.nxi;
    %             for i = 1:obj.Nss
    %                 %get the regulator equation solution                    
    %                 ind_pi = count + (1:nxn);
    %                 ind_th = count + nxn + (1:nxi);
    %                 Pi0{i} = sol0(ind_pi, :);
    %                 Th0{i} = sol0(ind_th, :);                    
    % 
    %                 count = count + nxn + nxi;
    % 
    %             end                                    
    %         end
    % 
    %         regulator_closed = struct('S', S, 'R', R, 'Pi', Pi0, ...
    %             'Th', Th0 );
    % 
    % 
    %     end
    % 
    %     %% use the model in synthesis
    % 
    %     function sys = get_model(obj, ind, vars_reg)
    %         %get_model
    %         %fetch the internal model (nominal) at mode 'ind'
    %         %
    %         %
    %         %with edits: allow for selection of model within feasible set
    % 
    %         %TODO: allow for parameterizations based on the variables
    % 
    % 
    %         if nargin < 3
    %             Phi = obj.Phi{ind};
    %             Gam = obj.Gam{ind};
    %         else
    %             Phi = vars_reg.Phi{ind};
    %             Gam = vars_reg.Gam{ind};
    %         end
    % 
    %         [nu, ns] = ssize(Gam);
    %         ny = ssize(Phi, 1);
    % 
    % 
    %         Am = obj.S;
    %         Bm = [zeros(ns, ny), eye(ns), zeros(ns, nu)];
    %         Cm = [-Gam; Phi];
    %         Dm = [zeros(nu, ny), zeros(nu, ns), eye(nu);
    %             eye(ny), zeros(ny, ns), zeros(ny, nu)];
    % 
    %         n = struct;
    %         n.nw = ny;
    %         n.nz = nu;
    %         n.nu = nu + ns;
    %         n.ny = ny;
    %         n.nzp = 0;
    %         n.nwp = 0;
    % 
    %         P = ss(Am, Bm, Cm, Dm, 1);
    %         sys = genplant(P, n);
    % 
    %     end
    % 
    %     function plant_model = connect_model(obj, plant, ind, rho)
    %         %connect the model (nominal regulator equation)
    % 
    %         if nargin < 4
    %             rho = 1;
    %         end
    % 
    %         if iscell(plant)
    %             plant_model = cell(obj.Nss, 1);
    %             for i = 1:obj.Nss
    %                 model = obj.get_model(i);
    %                 model_rho = rhotrafo(model, rho);
    %                 plant_model{i} = lft(plant{i}, model_rho);
    %             end
    %         else
    %             model = obj.get_model(ind);
    %             model_rho = rhotrafo(model, rho);
    %             plant_model = lft(plant, model_rho);
    %         end
    % 
    % 
    %     end
    % 
    % 
    % 
    end
end