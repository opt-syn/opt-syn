classdef regulator_lti < regulator_interface
    %REGULATOR_LTI Regulator for LTI systems
    

    
    methods
        function obj = regulator_lti(sys)
            %REGULATOR_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@regulator_interface(sys)
        end
        
%% Regulation check and internal models

        %solve the regulator equation for the network (Pi, Gamma, Phi)
        %   used for verification and internal model control
        %   find a parameterization of solutions of the regulator equation

        %construct the internal model

        %solve the regulator equation for the controller (Theta) given
        %(Gamma, Phi): this is verification

        function obj = form_internal_model(obj)
            %FORM_INTERNAL_MODEL create the internal model by solving the regulator
            %equation. Inputs are the system (P, bind, tracking, op)
            %
            %op is important for which oracles are equaltiy constarints and
            %which are inequality constraints

            %TODO: break this up into common routines


            %% break down the plant structure                      
            [A, B1, B2, C1, D11, D12, C2, D21, D22] = obj.sys.P.ss_zy_wu();

            
            Npre = obj.sys.get_consensus(obj.sys.op, obj.sys.bind);
            c = size(D11, 1)/length(obj.sys.bind); %coordinate lifts: change this later?
            N = kron(Npre, eye(c));
    
            [sN, dN] = size(N);
            n = obj.sys.P.nx;

            [Sbeta, Rbeta] = obj.sys.get_tracked_opt();

            if isempty(obj.sys.tracking)
                %constant tracking
                S = eye(size(N, 2)+c);
                R = S;


                reg_ans = [zeros(n, 1), -B1*N;  -ones(sN, 1), -D11*N];
                reg_mat = [A - eye(n), B2; C1, D12];


                null_basis = null(reg_mat, 'rational');
                try                    
                    sol0 = reg_mat \ reg_ans;
                catch
                    warning('Regulator equation cannot be solved')
                end
                nnull = size(null_basis, 2);

                Pi0 = sol0(1:n, :);
                Gam0 = sol0(n+1:end, :);        
                Phi0 = D21 * [zeros(sN, 1), N] + D22*Gam0 + C2*Pi0;

                if nnull
                    Pi_basis_pre = null_basis(1:n, :);
                    Gam_basis_pre = null_basis(n+1:end, :);
                    Phi_basis_pre = D22*Gam_basis_pre + C2*Pi_basis_pre;

                    Pi_basis = kron(Pi_basis_pre, eye(nnull));
                    Gam_basis = kron(Gam_basis_pre, eye(nnull));
                    Phi_basis = kron(Phi_basis_pre, eye(nnull));
                else
                    Pi_basis = [];
                    Gam_basis = [];
                    Phi_basis = [];
                end

            else
                %dynamical tracking, complicated Rbetalvester
                Sbeta = kron(Sbeta, eye(c));
                Rbeta = kron(Rbeta, eye(c));



                S = blkdiag(Sbeta, eye(dN));
                R = blkdiag(Rbeta, eye(dN));
                ns = size(S, 1);
                nr = size(Rbeta, 2);

                %use kronecker structure to explicitly solve the sylvester
                %equation
                reg_ans_mat = [zeros(n, size(Rbeta, 2)), -B1*N; -ones(sN, 1)*Rbeta, -D11*N];
    
                reg_ans = reshape(reg_ans_mat, [], 1);
                
                reg_mat_L = [A, B2; C1, D12];
                reg_mat_RR = S;
                reg_mat_RL = blkdiag(speye(size(A, 1)), sparse(sN, sN));

                reg_mat = kron(speye(ns), reg_mat_L) - kron(reg_mat_RR', reg_mat_RL);
                
                null_basis = null(reg_mat, 'rational');
                nnull = size(null_basis, 2);
                try                    
                    reg_sol_vec = reg_mat \ reg_ans;
                catch
                    warning('Regulator equation cannot be solved')
                end

                reg_sol = reshape(reg_sol_vec, size(reg_mat_L, 1), []);

                Pi0 = reg_sol(1:n, :);
                Gam0 = reg_sol(n+1:end, :);        
                Phi0 = D21 * [zeros(sN, nr), N] + D22*Gam0 + C2*Pi0;

                if nnull
                    %TODO: verify this 
                    Pi_basis = null_basis(1:n*ns, :);
                    Gam_basis = null_basis(ns*n+1:end, :);
                    Phi_basis = kron(D22, eye(nc))*Gam_basis + kron(C2, eye(nc))*Pi_basis;
    

                else
                    Pi_basis = [];
                    Gam_basis = [];
                    Phi_basis = [];
                end

            end


            obj.S = S;
            obj.R = R;
            obj.Pi = Pi0;
            obj.Gam = Gam0;
            obj.Phi = Phi0;
            obj.Pi_basis = Pi_basis;
            obj.Gam_basis = Gam_basis;
            obj.Phi_basis = Phi_basis;


        end


        function sys = get_model(obj, vars_reg)
            %get_model
            %fetch the internal model (nominal)
            %
            %
            %with edits: allow for selection of model within feasible set

            %TODO: allow for parameterizations based on the variables

            
            if nargin < 2
                Phi = obj.Phi;
                Gam = obj.Gam;
            else
                Phi = vars_reg.Phi;
                Gam = vars_reg.Gam;
            end

            [nu, ns] = ssize(Gam);
            ny = ssize(Phi, 1);


            Am = obj.S;
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

            P = ss(Am, Bm, Cm, Dm, 1);
            sys = genplant(P, n);

        end

        function plant_model = connect_model(obj, plant, rho)
            %connect the model (nominal regulator equation)

            if nargin < 3
                rho = 1;
            end

            model = obj.get_model();
            model_rho = rhotrafo(model, rho);
            plant_model = lft(plant, model_rho);
        end


        
        % function [regulator_closed] = verify_regulator(obj)
        %     %CHECK_REGULATOR does the closed-loop system obey the regulator
        %     %equation? (direct interconnection)
        % 
        %     regulator = obj.form_internal_model();
        % 
        %     [Ac, Bc, Cc, Dc] = ssdata(obj.K);
        % 
        %     %[Ac, Bc; Cc, Dc][Th; Gam] = [S; 0]Th'
        % 
        % 
        %                 reg_ans_mat = [zeros(n, size(Rbeta, 2)), -B1*N; -ones(sN, 1)*Rbeta, -D11*N];
        % 
        %         reg_ans = reshape(reg_ans_mat, [], 1);
        % 
        %         reg_mat_L = [A, B2; C1, D12];
        %         reg_mat_RR = S;
        %         reg_mat_RL = blkdiag(speye(size(A, 1)), sparse(sN, sN));
        % 
        % 
        % end

        function [regulator_closed] = check_regulator(obj)
            %CHECK_REGULATOR does the closed-loop system obey the regulator
            %equation? (direct interconnection)

            %TODO: simplify the computation based on the formed internal
            %model



            sys_cl = lft(obj.sys.P, obj.sys.K);

            A = sys_cl.A;
            B = sys_cl.Bw;
            C = sys_cl.Cz;
            D = sys_cl.Dzw;

            n = length(A);                

            
            Npre = obj.sys.get_consensus(obj.sys.op, obj.sys.bind);
            c = size(D, 1)/length(obj.sys.bind); %coordinate lifts: change this later?
            N = kron(Npre, eye(c));



            
            [sN, dN] = size(N);
            if isempty(obj.sys.tracking)
                %constant tracking
                
                S = eye(1+size(N, 2));
                R = S;

                reg_ans = [zeros(n, 1), -B*N;  -ones(sN, 1), -D*N];
                reg_mat = [A - eye(n); C];

                sol0 = reg_mat \ reg_ans;

                nx = obj.sys.P.nx;
                Pi0 = sol0(1:nx, :);
                Th0 = sol0(nx+1:end, :);        

               

            else

                %dynamical tracking, complicated Rbetalvester
                Sbeta = kron(obj.sys.tracking.Sbeta, eye(c));
                Rbeta = kron(obj.sys.tracking.Rbeta, eye(c));


                

                S = blkdiag(Sbeta, eye(dN));
                R = blkdiag(Rbeta, eye(dN));
                ns = size(S, 1);
                nr = size(Rbeta, 2);

                reg_ans_mat = [zeros(n, size(Rbeta, 2)), -B*N; -ones(sN, 1)*Rbeta, -D*N];
    
                reg_ans = reshape(reg_ans_mat, [], 1);
                
                reg_mat_L = [A; C];
                reg_mat_RR = S;
                reg_mat_RL = [speye(size(A, 1)); zeros(size(C))];


                reg_mat = kron(speye(ns), reg_mat_L) - kron(reg_mat_RR', reg_mat_RL);
                
                null_basis = null(reg_mat, 'rational');
                nnull = size(null_basis, 2);
                try                    
                    reg_sol_vec = reg_mat \ reg_ans;
                catch
                    warning('Regulator equation cannot be solved')
                end

                reg_sol = reshape(reg_sol_vec, [], ns);

                nx = obj.sys.P.nx;
                Pi0 = reg_sol(1:nx, :);
                Th0 = reg_sol(nx+1:end, :);        

            end


            
            % regulator_closed = regulator;
          regulator_closed = struct('S', S, 'R', R, 'Pi', Pi0, ...
              'Th', Th0 );

            

        end

        
    end
end

