classdef  opt_system
    %OPT_SYSTEM interconnection of network and operators
    %TODO: may be abstracted into an interface
    
    properties
        op; %a cell of operators (op_sim for simulation, op_? for analysis/synthesis)
        P;  %network
        K;  %controller
        bind; %which operators go to which output ports            
        tracking; %tracking of optimal solution (struct (S, R) by default)
                  %tracking of varying gradients requires LPV/periodic/switched
                  % methods, is a TODO
    end
    
    methods
        function obj = opt_system(op, P, K, bind, tracking)
            %OPT_SYSTEM constructor for the system
            obj.op = op;
            obj.P = P;
            obj.K = K;
            if nargin < 4
                s = length(obj.op);
                obj.bind = 1:s;
            else
                obj.bind = bind;
            end

            %assign identifiers to the operators
            for i = 1:length(op)
                obj.op{i}.id = i;
            end
            
            if nargin >= 5
                obj.tracking = tracking;
            end
            
        end    

        %% Dimension Counters
        function nss = Nss(obj)
            %NSS: number of subsystems
            nss = 1;
        end
       
        function dimn = n(obj)
            %n: number of states
            dimn = obj.nxn() + obj.nxi();
        end

        function dimn = nxn(obj)
            %nxn: number of states in network
            dimn = length(obj.P.A);
        end

        function dimn = nxi(obj)
            %nxi: number of states in controller
            dimn = length(obj.K.A);
        end
        

        %% getters
        function Pcurr = get_P(obj, param)
            %GET_P get the network P
            Pcurr = obj.P.ss();
        end
        
        function Kcurr = get_K(obj, param)
            %GET_K get the controller K
            Kcurr = obj.K;
        end

        function sys_alg = get_alg(obj, param)
            %close the loop of the algorithm
            if nargin < 2
                param = [];
            end
            Pcurr = obj.get_P(param);
            Kcurr = obj.get_K(param);
            sys_alg = lft(Pcurr, Kcurr);
        end

        function op_out = get_op(obj, i)
            %get the operator at index i
            op_out = obj.op{obj.bind(i)};
        end
        
        %% for simulation

        function [y, u] = get_internal_signals(obj, param, x_all, w_all)
            %extract the internal signals from the interconnection (y, u) 
            %using the well-posedness expression

            %Input:
            %   x_all:      all states of network and controller
            %   w_all:      all inputs to the network (except u)
            %   y_all:      all outputs to the network

            Kcurr = obj.get_K(param);
            Pcurr = obj.get_P(param);

            [nu, ny] = size(Kcurr.D);

            DK = Kcurr.D;
            DP = Pcurr.D((end-ny+1):end, (end-nu+1):end);


            nxi = length(Kcurr.A);
            nx = length(Pcurr.A);
            CyP = Pcurr.C((end-ny+1):end, :);            
            D21P = Pcurr.D((end-ny+1):end, 1:(end-nu));
            well_posed_mat = [eye(nu), -DP;
                              -DK, eye(ny)];
            

            nx = size(Pcurr.A, 1);
            nxi = size(Kcurr.A, 1);
            xN = x_all(1:(nx), :);
            xi = x_all((end-nxi+1):end, :);

            Cx = CyP*xN;
            Dw = D21P*w_all;
            Cxi = Kcurr.C * xi;

            sig_rhs = [Cx + Dw; Cxi];
            revert = well_posed_mat \ sig_rhs;


            y = revert(1:ny, :);
            u = revert((ny+1):end, :);
        end


        function mode_next = next_mode(obj, mode)
            %next mode in switching

            %TODO: is this actually used?
            mode_next = 1;
        end



        %% Regulation check and internal models

        %solve the regulator equation for the network (Pi, Gamma, Phi)
        %   used for verification and internal model control
        %   find a parameterization of solutions of the regulator equation

        %construct the internal model

        %solve the regulator equation for the controller (Theta) given
        %(Gamma, Phi): this is verification

        function regulator = form_internal_model(obj)
            %FORM_INTERNAL_MODEL create the internal model by solving the regulator
            %equation. Inputs are the system (P, bind, tracking, op)
            %
            %op is important for which oracles are equaltiy constarints and
            %which are inequality constraints


            %% break down the plant structure

            regulator = [];
            
            [A, B1, B2, C1, D11, D12, C2, D21, D22] = obj.P.ss_zy_wu();

            
            Npre = obj.get_consensus(obj.op, obj.bind);
            c = size(D11, 1)/length(obj.bind); %coordinate lifts: change this later?
            N = kron(Npre, eye(c));



            
            [sN, dN] = size(N);
            n = obj.P.nx;

            if isempty(obj.tracking)
                Sbeta = 1;
                Rbeta = 1;
            else
                Sbeta = obj.tracking.Sbeta;
                Rbeta = obj.tracking.Rbeta;
            end

            if isempty(obj.tracking)
                %constant tracking
                S = eye(size(N, 2)+1);
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
                    Pi_basis = null_basis(1:n, :);
                    Gam_basis = null_basis(n+1:end, :);
                    Phi_basis = D22*Gam_basis + C2*Pi_basis;
    

                else
                    Pi_basis = [];
                    Gam_basis = [];
                    Phi_basis = [];
                end

            end

            regulator = struct('S', S, 'R', R, 'Pi', Pi0, 'Gam', Gam0, 'Phi', Phi0, ...
                'Pi_basis', Pi_basis, 'Gam_basis', Gam_basis, 'Phi_basis',Phi_basis);


        end

        % 
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

            % regulator = obj.form_internal_model();


            sys_cl = lft(obj.P, obj.K);

            A = sys_cl.A;
            B = sys_cl.Bw;
            C = sys_cl.Cz;
            D = sys_cl.Dzw;

            n = length(A);                

            
            Npre = obj.get_consensus(obj.op, obj.bind);
            c = size(D, 1)/length(obj.bind); %coordinate lifts: change this later?
            N = kron(Npre, eye(c));



            
            [sN, dN] = size(N);
            if isempty(obj.tracking)
                %constant tracking
                
                S = eye(1+size(N, 2));
                R = S;

                reg_ans = [zeros(n, 1), -B*N;  -ones(sN, 1), -D*N];
                reg_mat = [A - eye(n); C];


                null_basis = null(reg_mat, 'rational');
                
                sol0 = reg_mat \ reg_ans;
                nnull = size(null_basis, 2);

                nx = obj.P.nx;
                Pi0 = sol0(1:nx, :);
                Th0 = sol0(nx+1:end, :);        

                if nnull
                    Pi_basis_pre = null_basis(1:nx, :);
                    Th_basis_pre = null_basis(nx+1:end, :);
    
                    Pi_basis = kron(Pi_basis_pre, eye(nnull));
                    Th_basis = kron(Th_basis_pre, eye(nnull));
                else
                    Pi_basis = [];
                    Th_basis = [];
                end

               

            else

                %dynamical tracking, complicated Rbetalvester
                Sbeta = kron(obj.tracking.Sbeta, eye(c));
                Rbeta = kron(obj.tracking.Rbeta, eye(c));


                

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

                nx = obj.P.nx;
                Pi0 = reg_sol(1:nx, :);
                Th0 = reg_sol(nx+1:end, :);        

                if nnull
                    Pi_basis_pre = null_basis(1:nx, :);
                    Th_basis_pre = null_basis(nx+1:end, :);
    
                    Pi_basis = kron(Pi_basis_pre, eye(nnull));
                    Th_basis = kron(Th_basis_pre, eye(nnull));
                else
                    Pi_basis = [];
                    Th_basis = [];
                end
            end


            
            % regulator_closed = regulator;
          regulator_closed = struct('S', S, 'R', R, 'Pi', Pi0, ...
              'Pi_basis', Pi_basis,...
              'Th', Th0,  'Th_basis', Th_basis);

            

        end

        function N = get_consensus(obj, op, bind)
            %GET_CONSENSUS create the consensus matrix
            %for the regulation condition

            %which operators are equality constraints
            % EQ = cellfun(@(e) e.EQUALITY, op);
            nop = length(op);
            EQ = zeros(1, nop, 'logical');
            for i = 1:nop
                EQ(i) = op{i}.EQUALITY;
            end

            
            if all(~EQ)
                s = length(op);
                N0 = [eye(s-1); -ones(1, s-1)];
                
            else  
                s = sum(~EQ);
                N0 = full(sparse(1:s, find(~EQ), ones(s, 1), nop, s));
                % N0 = [eye(s)];
            end

            %index based on the bind 
            nbind = length(bind);
            Bind = full(sparse(1:nbind, bind, ~EQ(bind), nbind, nop));


            N = Bind * N0;
        end        


    end
end

