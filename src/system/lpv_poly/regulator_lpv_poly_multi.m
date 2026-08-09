classdef regulator_lpv_poly_multi< regulator_interface
    %REGULATOR_LPV_POLY_MULTI Regulator for Polytopic LPV systems
    

    
    methods
        function obj = regulator_lpv_poly_multi(sys)
            %REGULATOR_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@regulator_interface(sys)

            % [obj.par, obj.pair] = obj.form_polytope();
        end
        % 
        % function [par, pair] = form_polytope(obj, sys)
        %     % Initialize parameters for the polytope representation
        %     par = obj.sys.par; % Placeholder for parameters
        %     pair = obj.sys.pair; % Placeholder for pairs
        % 
        % 
        % end

        

        %% dimension counters and indexers
        
        function NS = ns(obj)
            %NS: number of states of exosystem

            NS = length(obj.S{1});
        end

        function ns = Nss(obj)
            %NSS: number of parameters in the LPV model
            ns = size(obj.sys.par, 2);
        end

        function ns = Ncorner(obj)
            %NCORNER: number of corner/modes in the LPV model
            ns = size(obj.sys.par, 1);
        end

        function ns = Npair(obj)
            %NPAIR: number of transitions
            ns = size(obj.sys.pair, 1);
        end

        function P_accum = weight_sum(obj, P, parl)
            %weighted sum over the parameters

            P_accum = zeros(size(P{1}));
            for i = 1:length(parl)
                P_accum = P_accum + parl(i) * P{i};
            end
        end


        function P_corn = weight_sum_corners(obj)
            %weighted sum over each corner of the polytope
            P_corn = cell(obj.Ncorner, 1);
            for i = 1:obj.Ncorner   
                P_corn = obj.weight_sum(P, obj.par(i, :));
            end
        end       
    

        %% routines to compute the internal model
        function [reg_mat_all, reg_ans_all] = reg_sys_all(obj)
            %assemble the regulator equation system

            reg_mat_all = [];
            reg_ans_all = [];

            I = eye(obj.Nss);
            for i = 1:obj.Npair
    
                %form the transition constraint
                pair_curr = obj.sys.pair(i, :);
                par_curr = pair_curr(:, 1:obj.Nss);
                par_next = pair_curr(:, (obj.Nss+1):end);



                reg_mat_S = [];
                reg_mat_dyn = [];
                reg_ans = 0;

                [reg_mat_dyn_curr, reg_ans_curr] = obj.reg_sys_indiv(par_curr);    
                reg_mat_S_next = obj.reg_sys_next(par_next);


                %take kroneckers to find the per-subsystem contributions
                reg_ans = reg_ans_curr;
                reg_mat_dyn = kron(par_curr, reg_mat_dyn_curr);
                reg_mat_S = kron(par_next, reg_mat_S_next);
                

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
                parl = I(i, :);

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
                parl = I(i, :);
                null_curr = squeeze(null_basis(count + (1:(nxn + nxu)), :, :));
                [Pi_basis{i}, Gam_basis{i}, Phi_basis{i}] = null_reg_index(obj, null_curr, parl);
                count = count + nxn + nu;
            end
        end

        %% acquire the models
        function sys = get_model(obj, parl, vars_reg)
            %get_model
            %fetch the internal model (nominal)
            %
            %parl: specific parameter value
            %
            %
            %these internal models are affine in theta
            %
            %maybe they should also be scheduled?           
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

            S = obj.weight_sum(obj.S, parl);
            Phi = obj.weight_sum(Phi, parl);
            Gam = obj.weight_sum(Gam, parl);

            sys =obj.fetch_model(S, Phi, Gam);
            % sys = cell(obj.Nss, 1);
            % for i = 1:obj.Nss
            %     sys{i} =obj.fetch_model(obj.S{i}, Phi{i}, Gam{i});
            % end
            

        end

        function plant_model = connect_model(obj, plant, ind_par, rho)
            %connect the model (nominal regulator equation)

            if nargin < 3
                rho = 1;
            end

            if ~iscell(plant)
                plant0 = plant;
                plant = cell(obj.Nss, 1);
                for i = 1:obj.Nss
                    plant{i} = plant0;
                end
            end


            %form the model interconnection
            I = eye(obj.Nss);
            % for i = 1:obj.Nss
                parl = I(ind_par, :);
                model = obj.get_model(parl);
                model_rho = rhotrafo(model, rho);
                plant_model = lft(plant, model_rho);
            % end

        end
        
    end
end

