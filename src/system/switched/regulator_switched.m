classdef regulator_switched < regulator_interface
    %REGULATOR_SWITCHED Regulator for switched systems 
    %
    

    % [x(k+1)] = [A(mode(k))    Bd(mode(k))    Bu(mode(k))  ][x(k)]   state transition
    % [e(k)  ] = [Ce(mode(k))   Ded(mode(k))   Deu(mode(k)) ][d(k)]   output to  regulated error    
    % [u(k) ] =  [Cy(mode(k))   Dyd(mode(k))   Dyu(mode(k)) ][u(k)]   output to controller    


    %The system obeys a switching logic: restricted mode transitions 
    %from mode(k) to mode(k+1) based on a switching graph (adjacency matrix
    % sys.adj)
    %
    %instances of these algorithms include optimization algorithms under
    %a-priori-unknown time-varying delays

    methods
        function obj = regulator_switched(sys)
            %REGULATOR_SWITCHED build the regulator
            
            obj@regulator_interface(sys)
        end

        function ns = Nss(obj)
            %NSS: number of subsystems
            ns = obj.sys.Nss;
        end

        %% form the internal model
        function [reg_mat_all, reg_ans_all] = reg_sys_all(obj)
            %assemble the regulator equation system

            reg_mat_all = [];
            reg_ans_all = [];

            %alignment across all arcs
            [src, dst] = obj.sys.get_arcs();
            Narcs = length(src);

            

            for i = 1:Narcs

                %form the transition constraint
                % pair_curr = obj.sys.pair(i, :);
                % par_curr = pair_curr(:, 1:obj.Nss);
                % par_next = pair_curr(:, (obj.Nss+1):end);

                par_curr = struct('mode', src(i));

                [reg_mat_dyn_L_curr, reg_mat_dyn_R_curr, reg_ans_curr] = obj.reg_sys_indiv(par_curr);    
                [reg_mat_S_L_next, reg_mat_S_R_next] = obj.reg_sys_next(par_curr);
                
                %take kroneckers to find the per-subsystem contributions
                reg_ans = reg_ans_curr;
                pcurr = full(sparse(1, src(i), 1, 1, obj.Nss));
                pnext = full(sparse(1, dst(i), 1, 1, obj.Nss));

                reg_mat_dyn_L = kron(pcurr, reg_mat_dyn_L_curr);
                reg_mat_S_L = kron(pnext, reg_mat_S_L_next);


                reg_mat_dyn = kron(reg_mat_dyn_R_curr', reg_mat_dyn_L);
                reg_mat_S = kron(reg_mat_S_R_next', reg_mat_S_L);

                %append the answer
                reg_ans_all = [reg_ans_all; reg_ans];

                %load in to the matrix               
                % reg_mat_expand = ;

                reg_mat_all = [reg_mat_all; reg_mat_dyn - reg_mat_S];
            end
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
                        %the first term is constant. Then the rest are
                        %deviations.
                        if i > 1
                            Isupp = zeros(dN);
                        else
                            Isupp = eye(dN);
                        end
                        S{i} = blkdiag(Sbeta{i}, Isupp);
                        R{i} = blkdiag(Rbeta{i}, Isupp);
                    end
                else
                    S = blkdiag(Sbeta, eye(dN));
                    R = blkdiag(Rbeta, eye(dN));
                end
            end
        end

        function NS = ns(obj)
            %NS number of states of exosystem

            if iscell(obj.S)
                NS = length(obj.S{1});
            else
                NS = length(obj.S);
            end
        end



        function [Pi, Gam, Phi] = sol_reg_all(obj, reg_sol)
            %recover the solution to the regulator equation system

            ns = obj.ns;
            reg_sol = reshape(reg_sol, [], ns);

            nxn = obj.sys.nxn;
            nu = obj.sys.nu;

            count = 0;
            I = eye(obj.Nss);
            Pi = cell(obj.Nss, 1);
            Gam = cell(obj.Nss, 1);
            Phi = cell(obj.Nss, 1);
            for i = 1:obj.Nss

                reg_sol_curr = reg_sol(count + (1:(nxn + nu)), :);
                % parl = I(i, :);
                parl = struct('mode', i);

                count = count + nxn + nu;
                [Pi{i}, Gam{i}, Phi{i}] = obj.sol_reg_index(reg_sol_curr, parl);
            end

        end

        function [Pi_basis, Gam_basis, Phi_basis] = null_reg(obj, null_basis);
            %NULL_REG a nullspace indexer (altogether)


            if nargin < 3
                param =[];
            end

            Pi_basis = cell(obj.Nss, 1);
            Gam_basis = cell(obj.Nss, 1);
            Phi_basis = cell(obj.Nss, 1);
            I = eye(obj.Nss);

            count = 0;
            nxn = obj.sys.nxn;
            nu = obj.sys.nu;

            for i = 1:obj.Nss
                % parl = I(i, :);
                parl = struct('mode', i);

                null_curr = squeeze(null_basis(count + (1:(nxn + nu)), :, :));
                [Pi_basis{i}, Gam_basis{i}, Phi_basis{i}] = null_reg_index(obj, null_curr, parl);
                count = count + nxn + nu;
            end
        end


        %% check the regulator equations

        function [reg_mat_all, reg_ans_all] = reg_K_sys_all(obj)
            %assemble the closed-loop regulator equation system 
            %
            %Returns: 
            %   reg_mat:    matrix for regulator equation
            %   reg_ans:    vector for regulator equation solution
           

            reg_mat_all = [];
            reg_ans_all = [];

            %alignment across all arcs
            [src, dst] = obj.sys.get_arcs();
            Narcs = length(src);

            for i = 1:Narcs

                %form the transition constraint  
                par_curr = struct('mode', src(i));
                par_next = struct('mode', dst(i));

                [reg_mat_dyn_L_curr, reg_mat_dyn_R_curr, reg_ans_curr] = obj.reg_K_sys_indiv(par_curr);    
                [reg_mat_S_L_next, reg_mat_S_R_next] = obj.reg_K_sys_next(par_curr);

                %take kroneckers to find the per-subsystem contributions
                reg_ans = reg_ans_curr;
                pcurr = full(sparse(1, src(i), 1, 1, obj.Nss));
                pnext = full(sparse(1, dst(i), 1, 1, obj.Nss));

                reg_mat_dyn_L = kron(pcurr, reg_mat_dyn_L_curr);
                reg_mat_S_L = kron(pnext, reg_mat_S_L_next);


                reg_mat_dyn = kron(reg_mat_dyn_R_curr', reg_mat_dyn_L);
                reg_mat_S = kron(reg_mat_S_R_next', reg_mat_S_L);

 

                %append the answer
                reg_ans_all = [reg_ans_all; reg_ans];

                %load in to the matrix               
                % reg_mat_expand = ;

                reg_mat_all = [reg_mat_all; reg_mat_dyn - reg_mat_S];
            end
        end

        function [reg_mat_dyn_L, reg_mat_dyn_R, reg_ans] = reg_K_sys_indiv(obj, param)
            %control regulator equation checks (closed-loop)

            if nargin < 2
                param = [];
            end
            Kcurr = obj.get_K(param);

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
            reg_mat_dyn_L = [Ak; Ck];

            reg_mat_dyn_R = speye(ns);

            if ~isempty(obj.Gam_basis)
                nnull = size(obj.Gam_basis, 3);

                Gam_contract = tensorprod([Bk; Dk], obj.Gam_basis, 2, 1);
                Phi_contract = tensorprod([Bk; Dk], obj.Phi_basis, 2, 1);
            end

        end

        function Pi_out = get_Pi(obj, param)
            if nargin == 2
                Pi_out = obj.Pi{param.mode};
            else
                Pi_out = obj.Pi;
            end
        end

        function Gam_out = get_Gam(obj, param)            
            if nargin == 2
                Gam_out = obj.Gam{param.mode};
            else
                Gam_out = obj.Gam;
            end
        end

        function Phi_out = get_Phi(obj, param)
            if nargin == 2
                Phi_out = obj.Phi{param.mode};
            else
                Phi_out = obj.Phi;
            end
        end

        function [Pi0, Gam0, Phi0, Th0] = sol_K_reg_all(obj, reg_sol)
            %recover the solution to the regulator equation system

            ns = obj.ns;
            reg_sol = reshape(reg_sol, [], ns);


            ns = obj.ns;
            reg_sol = reshape(reg_sol, [], ns);

            nxi = obj.sys.nxi;
            

            count = 0;
            
            Pi0 = cell(obj.Nss, 1);
            Gam0 = cell(obj.Nss, 1);
            Phi0 = cell(obj.Nss, 1);
            Th0 = cell(obj.Nss, 1);

            for i = 1:obj.Nss

                reg_sol_curr = reg_sol(count + (1:nxi), :);
                % parl = I(i, :);
                parl = struct('mode', i);

                count = count + nxi;
                [Pi0{i}, Gam0{i}, Phi0{i}, Th0{i}] = sol_K_reg_index(obj, reg_sol_curr, parl);

            end

        end

    

        %% use the model in synthesis

        function sys = get_model(obj, ind, vars_reg)
            %get_model
            %fetch the internal model (nominal) at mode 'ind'
            %with edits, allow for selection of model within feasible set
            %Args:
            %   vars_reg:   variables of the problem        (regulator)            
            %Return:
            %   model: the full-order internal model
            

            %TODO: allow for parameterizations based on the variables

            
            if nargin < 3
                Phi = obj.Phi{ind};
                Gam = obj.Gam{ind};
            else
                Phi = vars_reg.Phi{ind};
                Gam = vars_reg.Gam{ind};
            end

            [nu, ns] = ssize(Gam);
            ny = ssize(Phi, 1);


            if iscell(obj.S)
                Am = obj.S{ind};
            else
                Am = obj.S;
            end
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

        function plant_model = connect_model(obj, plant, ind, rho)
            %connect the model (nominal regulator equation)
            %
            %Args:
            %   plant: original system
            %   ind:    index to examine
            %   rho:    exponential weighting
            %Return:
            %   plant_model: plant and model together

            if nargin < 4
                rho = 1;
            end

            if iscell(plant)
                plant_model = cell(obj.Nss, 1);
                for i = 1:obj.Nss
                    model = obj.get_model(i);
                    model_rho = rhotrafo(model, rho);
                    plant_model{i} = lft(plant{i}, model_rho);
                end
            else
                model = obj.get_model(ind);
                model_rho = rhotrafo(model, rho);
                plant_model = lft(plant, model_rho);
            end

            
        end

        function vars_reg = create_vars(obj)
            %CREATE_VARS: create variables that parameterize the nullspace
            %Args:
            %   param_null(bool): should the nullspace be searched (as
            %   variables)
            %
            %Returns:
            %   vars_reg: structure with fields (Pi, Gam, Phi)
            
            vars_reg.Pi = obj.Pi;
            vars_reg.Gam = obj.Gam;
            vars_reg.Phi = obj.Phi;

            

            
            if isempty(obj.Gam_basis) || isempty(obj.Gam_basis{1})
                vars_reg.eta= [];
            % else
            %     nbasis = size(obj.Gam_basis{1}, 3);                
            %     eta =  lmim('reg_param', nbasis, 1, 'full');
            % 
            % 
            %     vars_reg.Pi = cell(obj.Nss, 1);
            %     vars_reg.Gam = cell(obj.Nss, 1);
            %     vars_reg.Phi = cell(obj.Nss, 1);
            %     for i = 1:nbasis
            %         for j = 1:obj.Nss
            %             eta_curr = lmim_index(eta, i, 1);
            % 
            %             vars_reg.Pi{j} = vars_reg.Pi{j} + obj.Pi_basis{j}(:, :, i) * eta_curr;
            %             vars_reg.Gam{j} = vars_reg.Gam{j} + obj.Gam_basis{j}(:, :, i) * eta_curr;
            %             vars_reg.Phi{j} = vars_reg.Phi{j} + obj.Phi_basis{j}(:, :, i) * eta_curr;
            %         end
            %     end
            % 
            %     vars_reg.eta = eta;
            % 
                
            end
        end
        
    




    end
end