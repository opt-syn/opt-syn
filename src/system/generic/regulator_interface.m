classdef regulator_interface
    %REGULATOR_INTERFACE 
    %
    % A regulator employed for the synthesis of optimization algorithms
    %
    % This regulator carries the necessary internal model required for
    % convergence of optimization algorithms (constant shift of optimal 
    % solution), as well as extra states associated with position of the
    % unknown disturbance. This regulator is adjoined to the network
    % dynamics in synthesis of algorithms.
    %
    % The regulator allows for the perfect tracking of optimal solutions
    % (if known) by the tracking dynamics (Sbeta, Rbeta) in (sys).
    %
    %
    % The regulator is specialized for a specific kind of system
    %
    % 
    %
    %Open issues:
    % multiple solutions to the regulator equations
    
    properties
        sys;

        %the exosystem for the optimal solution and operators
        S;
        R;

        %nominal solution to the regulator equations
        Pi;
        Gam;
        Phi;

        %freedom in solving the regulator equations
        Pi_basis;
        Gam_basis;
        Phi_basis;
    end
    
    methods
        function obj = regulator_interface(sys)
            %REGULATOR_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.sys = sys;

            obj = obj.form_internal_model();
        end

        function NS = ns(obj)
            %NS: number of states of exosystem

            NS = length(obj.S);
        end

        function sys = fetch_model(obj, S, Phi, Gam)
            %FETCH_MODEL: fetch an internal model from data
            [nu, ns] = ssize(Gam);
            ny = ssize(Phi, 1);


            Am = S;
            Bm = [zeros(ns, ny), eye(ns), zeros(ns, nu)];
            Cm = [-Gam; Phi];
            Dm = [zeros(nu, ny), zeros(nu, ns), eye(nu);
                eye(ny), zeros(ny, ns), zeros(ny, nu)];

            n = struct;
            n.nw = ny;
            n.nz = nu;
            n.nu = nu + ns;
            n.ny = ny;
            n.nzp = 0;
            n.nwp = 0;

            if isnumeric(Phi) && isnumeric(Gam) && isnumeric(S)
                P = ss(Am, Bm, Cm, Dm, 1);

                statename = {};
                inputname = {};
                outputname = {};
                for i = 1:ns
                    statename{i} = ['v', num2str(i)];
                end

                for i = 1:ny+ns+nu
                    if i <= ny
                        inputname{i} = ['y', num2str(i)];
                    elseif i <= ns + ny
                        inputname{i} = ['um1_', num2str(i - ns)];
                    else
                        inputname{i} = ['um2_', num2str(i-ns-ny)];
                    end
                end

                for i = 1:ny+nu
                    if i <= nu
                        outputname{i} = ['u', num2str(i)];
                    else
                        outputname{i} = ['ym', num2str(i - nu)];
                    end
                end
                   
                P.StateName = statename;
                P.InputName = inputname;
                P.OutputName = outputname;
                % P.InputNames = 
            else
                P = sdpss(Am, Bm, Cm, Dm);
            end

            sys = genplant(P, n);


        end

        
        %% solve the regulator equations (open loop)
        function obj = form_internal_model(obj)
            %FORM_INTERNAL_MODEL create the internal model by solving the regulator
            %equation. Inputs are the system (P, bind, tracking, op)
            %
            %op is important for which oracles are equaltiy constarints and
            %which are inequality constraints

            %TODO: break this up into common routines
 
                
            [S, R] = obj.exosystem();
            obj.S = S;
            obj.R = R;

            [reg_mat, reg_ans] = obj.reg_sys_all();
            
            
            reg_sol= lsqminnorm(reg_mat, reg_ans);
            % catch
            reg_err = reg_mat * reg_sol - reg_ans;
            if norm(reg_err) > 1e-8
            
                warning('Regulator equation cannot be solved')
            end

                
            [Pi0, Gam0, Phi0] = obj.sol_reg_all(reg_sol);
            [Pi_basis, Gam_basis, Phi_basis] = obj.null_reg_all(reg_mat);

            
            obj.Pi = Pi0;
            obj.Gam = Gam0;
            obj.Phi = Phi0;
            obj.Pi_basis = Pi_basis;
            obj.Gam_basis = Gam_basis;
            obj.Phi_basis = Phi_basis;       
        end

        
        function [Pi, Gam, Phi] = sol_reg_all(obj, reg_sol)
            %recover the solution to the regulator equation system

            ns = obj.ns;
            reg_sol = reshape(reg_sol, [], ns);
            [Pi, Gam, Phi] = obj.sol_reg_index(reg_sol);
        end

        % individual regulator terms
        function N = get_consensus(obj)
            %get the consensus matrix
            Npre = obj.sys.get_consensus(obj.sys.op, obj.sys.bind);
            c = obj.sys.op{1}.c; %coordinate lifts: change this later?
            N = kron(Npre, eye(c));
        end

        function[Pi, Gam, Phi] = sol_reg_index(obj, reg_sol, param)
            %index the solution to the regulator equation
            if nargin < 3
                param =[];
            end

            
            
            n = obj.sys.nxn;
            Pi = reg_sol(1:n, :);
            Gam = reg_sol(n+1:end, :);        
            Phi = obj.compute_Phi(Pi, Gam);
        end

        function [Phi] = compute_Phi(obj, Pi, Gam, param)
            if nargin < 4
                param = [];
            end
            [A, B1, B2, C1, D11, D12, C2, D21, D22] = obj.sys.ss_zy_wu(param);
            [Bd, Ded, Dyd] = obj.d_influence(param);
            Phi = Dyd + D22*Gam + C2*Pi;
        end


        function [Pi_basis, Gam_basis, Phi_basis] = null_reg_all(obj, reg_mat)
            %NULL_REG nullspace of the regulator equations: freedom to move            
            
            null_basis = null(reg_mat, 'rational');
            nnull = size(null_basis, 2);

            if nnull      
                ns = obj.ns;
  
                null_basis = reshape(null_basis, [], ns, nnull);
                [Pi_basis, Gam_basis, Phi_basis] = obj.null_reg(null_basis);
            else
                Pi_basis = [];
                Gam_basis = [];
                Phi_basis = [];
            end
                

        end
        
        function [Pi_basis, Gam_basis, Phi_basis] = null_reg(obj, null_basis);
            %NULL_REG a nullspace indexer (altogether)


            if nargin < 3
                param =[];
            end

            [Pi_basis, Gam_basis, Phi_basis] = null_reg_index(obj, squeeze(null_basis), param);
        end

        function [Pi_basis, Gam_basis, Phi_basis] = null_reg_index(obj, null_basis, param);
            %NULL_REG_INDEX a  nullspace indexer for each subsystem
                       
            ns = obj.ns;
            n = obj.sys.nxn;

            [A, B1, B2, C1, D11, D12, C2, D21, D22] = obj.sys.ss_zy_wu(param);

            Pi_basis =null_basis(1:n, :, :);
            Gam_basis = null_basis((n+1):end, :, :);
            Phi_basis = tensorprod(D22, Gam_basis, 2, 1) + tensorprod(C2, Pi_basis, 2, 1);   

        end

        function [reg_mat, reg_ans] = reg_sys_all(obj)
            %assemble the regulator equation system


            [reg_mat_dyn, reg_ans] = obj.reg_sys_indiv();            

            reg_mat_S = obj.reg_sys_next();
            
            reg_mat = reg_mat_dyn - reg_mat_S;
        end

        function reg_mat_S = reg_sys_next(obj, param)
            %REG_SYS the system to be regulated, forming the regulator
            %equations
            %the system to be regulated (the next state/sylvester expression)

            if nargin < 2

                param = [];
            end

            [S, R] = obj.exosystem(param);
            N = obj.get_consensus();
            [sN, dN] = size(N);

            reg_mat_RR = S;
            n = obj.sys.nxn;
            nu = obj.sys.nu;
            reg_mat_RL = [speye(n), sparse(n, nu);
                         sparse(sN, n), sparse(sN, nu)];

            
            reg_mat_S = kron(reg_mat_RR', reg_mat_RL);

        end

        function [reg_mat_dyn, reg_ans] = reg_sys_indiv(obj, param)
            %REG_SYS the system to be regulated, forming the regulator
            %equations
            %the system to be regulated (dynamics expression)
            if nargin < 2

                param = [];
            end

            
            

            [Sbeta, Rbeta] = obj.get_tracked_opt(param);
            [A, B1, B2, C1, D11, D12, C2, D21, D22] = obj.sys.ss_zy_wu(param);

            [Bd, Ded, Dyd] = obj.d_influence(param);
            [S, R] = obj.exosystem(param);
            ns = length(S);

            N = obj.get_consensus();
            [sN, dN] = size(N);
            %form the system
    
            reg_ans_mat = -[Bd; Ded];
    
            reg_ans = reshape(reg_ans_mat, [], 1);
            
            reg_mat_L = [A, B2; C1, D12];
            reg_mat_dyn = kron(speye(ns), reg_mat_L);


        end

        % 
        function [Bd, Ded, Dyd] = d_influence(obj, param)
           %D_INFLUENCE how does the system get affected by the
           %disturbance?
           %
           %
           %used for the internal model computation
           if nargin < 2
               param =[];
           end

            N = obj.get_consensus();
            [sN, dN] = size(N);
            n = obj.sys.P.nx;
            c = obj.sys.op{1}.c;

            [Sbeta, Rbeta] = obj.get_tracked_opt(param);
            [A, B1, B2, C1, D11, D12, C2, D21, D22] = obj.sys.ss_zy_wu(param);

            Bd = [zeros(n, size(Rbeta, 2)), B1*N];
            Ded = [kron(ones(sN/c, 1), eye(c))*Rbeta, D11*N];
            Dyd = [zeros(size(D21, 1), c)*Rbeta, D21*N];

        end

        function saug = sys_regulated_aug(obj)
            %the enriched system with the disturbance channel
            saug = [];
        end

        %% fetch the exosystem
        function [S, R] = exosystem(obj, param)

            %get the exosystem at each mode/internal model
            N = obj.get_consensus();
            [sN, dN] = size(N);


            if nargin == 2
                [Sbeta, Rbeta] = obj.sys.get_tracked_opt(param);                   

                S = blkdiag(Sbeta, eye(dN));
                R = blkdiag(Rbeta, eye(dN));
            else
                [Sbeta, Rbeta] = obj.sys.get_tracked_opt();
                if iscell(Sbeta)
                    S = cell(size(Sbeta));
                    R = cell(size(Rbeta));
                    for i = 1:numel(S)
                        S{i} = blkdiag(Sbeta{i}, eye(dN));
                        R{i} = blkdiag(Rbeta{i}, eye(dN));
                    end
                else
                    S = blkdiag(Sbeta, eye(dN));
                    R = blkdiag(Rbeta, eye(dN));
                end
            end
        end

        function [Sbeta, Rbeta] = get_tracked_opt(obj, param)

            %get the part of the exosystem corresponding to tracking the
            %optimal solution at each mode/internal model            

            if nargin == 2
                [Sbeta, Rbeta] = obj.sys.get_tracked_opt(param);                   

            else
                [Sbeta, Rbeta] = obj.sys.get_tracked_opt();
            end
        end

        %% check the regulator equations (closed-loop)

        function [reg_cl] = check_regulator(obj)
            %check the regulator equation for a specific system
            
            N = obj.get_consensus();

            

            [reg_mat, reg_ans] = obj.reg_K_sys_all();
            
            % try                    
                reg_sol= lsqminnorm(reg_mat, reg_ans);
            % catch
            reg_err = reg_mat * reg_sol - reg_ans;
            if norm(reg_err) > 1e-8
                cerror('Regulator equation cannot be solved')
            end
            % end

            [Pi, Gam, Phi, Th] = obj.sol_K_reg_all(reg_sol);

            reg_cl = struct;
            reg_cl.S = obj.S;
            reg_cl.R = obj.R;
            reg_cl.Pi = Pi;
            reg_cl.Gam = Gam;
            reg_cl.Phi = Phi;
            reg_cl.Th = Th;

        end

        function Pi_out = get_Pi(obj, param)
            Pi_out = obj.Pi;
        end

        function Gam_out = get_Gam(obj, param)            
            Gam_out = obj.Gam;
        end

        function Phi_out = get_Phi(obj, param)
            Phi_out = obj.Phi;
        end

        function [reg_mat, reg_ans] = reg_K_sys_all(obj)
            %assemble the regulator equation system

            [reg_mat_dyn, reg_ans] = obj.reg_K_sys_indiv();            

            reg_mat_S = obj.reg_K_sys_next();
            reg_mat = reg_mat_dyn - reg_mat_S;
        end

        function [reg_mat_S] = reg_K_sys_next(obj, param)
            %control regulator equation checks
            %next expression
            if nargin < 2
                param = [];
            end

            [S, R] = obj.exosystem(param);

            N = obj.get_consensus();
            [sN, dN] = size(N);

            reg_mat_RR = S;
            nxi = obj.sys.nxi;
            reg_mat_RL = [speye(nxi); sparse(sN, nxi)];


            reg_mat_S = kron(reg_mat_RR', reg_mat_RL);

        end

        function [reg_mat_dyn, reg_ans] = reg_K_sys_indiv(obj, param)
            %control regulator equation checks
            %dynamics expression
            %this excludes nullspace (for now)

            if nargin < 2
                param = [];
            end
            Kcurr = obj.sys.get_K(param);

            [Ak, Bk, Ck, Dk] = ssdata(Kcurr);
            [S, R] = obj.exosystem(param);

            N = obj.get_consensus();
            [sN, dN] = size(N);
    

            %answer to regulator equation
            Gam = obj.get_Gam(param);
            Phi= obj.get_Phi(param);

            reg_ans = reshape([-Bk*Phi; Gam - Dk*Phi], [], 1);

            
            %system for 
            nxi = size(Ak, 1);
            ns = obj.ns;
            % nu = 
            reg_mat_L = [Ak; Ck];

            reg_mat_dyn = kron(speye(ns), reg_mat_L);   

            if ~isempty(obj.Gam_basis)
                nnull = size(obj.Gam_basis, 3);
                
                Gam_contract = tensorprod([Bk; Dk], obj.Gam_basis, 2, 1);
                Phi_contract = tensorprod([Bk; Dk], obj.Phi_basis, 2, 1);
            end

        end

        function [Pi0, Gam0, Phi0, Th0] = sol_K_reg_all(obj, reg_sol)
            %recover the solution to the regulator equation system

            ns = obj.ns;
            reg_sol = reshape(reg_sol, [], ns);
            [Pi0, Gam0, Phi0, Th0] = obj.sol_K_reg_index(reg_sol, []);
        end

        function[Pi, Gam, Phi, Th] = sol_K_reg_index(obj, reg_sol, param)
            %index the solution to the regulator equation
            if nargin < 3
                param =[];
            end

            Pi  = obj.get_Pi(param);
            Gam = obj.get_Gam(param);
            Phi = obj.get_Phi(param);

            ns = obj.ns;
            nxi = obj.sys.nxi;
           

            %todo: finish this, including nullspace entries
            Th = reshape(reg_sol(1:(ns*nxi)), nxi, ns);
            
            
        end



        %% incorporate into optimization
        function vars_reg = create_vars(obj)
            %CREATE_VARS: create variables that parameterize the nullspace
            vars_reg.Pi = obj.Pi;
            vars_reg.Gam = obj.Gam;
            vars_reg.Phi = obj.Phi;

            

            
%             if isempty(obj.Gam_basis)==0
            if false
                nbasis = size(obj.Gam_basis, 3);                
                eta = lmim('reg_param', nbasis, 1, 'full');


                for i = 1:nbasis
                    eta_curr = lmim_index(eta, i, 1);

                    vars_reg.Pi = vars_reg.Pi + obj.Pi_basis(:, :, i) * eta_curr;
                    vars_reg.Gam = vars_reg.Gam + obj.Gam_basis(:, :, i) * eta_curr;
                    vars_reg.Phi = vars_reg.Phi + obj.Phi_basis(:, :, i) * eta_curr;
                end

                vars_reg.eta = eta;
            else
                vars_reg.eta= [];
            end
        end
        
    end

    methods(Abstract)
        % form_internal_model(obj)
        % check_regulator(obj)
        % get_model(obj)
    end
end

