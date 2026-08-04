classdef regulator_lpv_poly< regulator_interface
    %REGULATOR_LPV_POLY Regulator for Polytopic LPV systems
    %
    %
    %This is a constant solution to the regulator equation
    %the regulator equation conditions can get bilinear and tricky if 
    %parameter-dependent solutions are allowed. Mode-dependent solutions
    %are fine for switched systems though.

    %[A(th),  Bd(th),  Bu(th) ][Pi]  = [Pi S(th_next)]
    %[Ce(th), Ded(th), Deu(th)][I]   = [0]
    %[Cy(th), Dyd(th), Dyu(th)][Gam] = [Phi(th)]
    %
    %
    %[Pi, Gam] are independent of theta (makes the formulation easier to
    %solve, but does add conservatism in optimization algorithm design)


    properties
        S_corner = [];
    end
    methods
        function obj = regulator_lpv_poly(sys)
            %REGULATOR_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@regulator_interface(sys)    


            obj.S_corner = obj.weight_sum_corners(obj.S);
        end

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


        function P_corn = weight_sum_corners(obj, P)
            %weighted sum over each corner of the polytope
            P_corn = cell(obj.Ncorner, 1);
            for i = 1:obj.Ncorner   
                P_corn{i} = obj.weight_sum(P, obj.sys.par(i, :));
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


                reg_ans = reg_ans_curr;

                %append the answer
                reg_ans_all = [reg_ans_all; reg_ans];

                %load in to the matrix               
                % reg_mat_expand = ;
                
                %parameter-independent regulator equation solution
                reg_mat_all = [reg_mat_all; reg_mat_dyn_curr - reg_mat_S_next];

                % reg_mat_all = [reg_mat_all; reg_mat_dyn - reg_mat_S];
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
            Pi = cell(obj.Ncorner, 1);
            Gam = cell(obj.Ncorner, 1);
            Phi = cell(obj.Ncorner, 1);

            % Pi0 = reg_sol(1:nxn, :);
            % Gam0 = reg_sol(nxn + (1:nu), :);
            reg_sol_curr = reg_sol(count + (1:(nxn + nu)), :);
            for i = 1:obj.Ncorner                                              
                % parl = I(i, :);
                % count = count + nxn + nu;
                [Pi{i}, Gam{i}, Phi{i}] = obj.sol_reg_index(reg_sol_curr, obj.sys.par(i, :));
            end

        end

        function [reg_mat_all, reg_ans_all] = reg_K_sys_all(obj)
            %assemble the regulator equation system

            reg_mat_all = [];
            reg_ans_all = [];

            
            for i = 1:obj.Npair
    
                %form the transition constraint
                pair_curr = obj.sys.pair(i, :);
                par_curr = pair_curr(:, 1:obj.Nss);
                par_next = pair_curr(:, (obj.Nss+1):end);



                reg_mat_S = [];
                reg_mat_dyn = [];
                reg_ans = 0;

                [reg_mat_dyn_curr, reg_ans_curr] = obj.reg_K_sys_indiv(par_curr);    
                reg_mat_S_next = obj.reg_sys_next(par_next);


                reg_ans = reg_ans_curr;

                %append the answer
                reg_ans_all = [reg_ans_all; reg_ans];

                %load in to the matrix               
                % reg_mat_expand = ;
                
                %parameter-independent regulator equation solution
                reg_mat_all = [reg_mat_all; reg_mat_dyn_curr - reg_mat_S_next];

                % reg_mat_all = [reg_mat_all; reg_mat_dyn - reg_mat_S];
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
        function sys = get_model(obj, parl)
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
            
            % if nargin < 3
            ind = find(parl);
                Phi = obj.Phi{ind};
                Gam = obj.Gam{ind};
            % else
                % Phi = vars_reg.Phi;
                % Gam = vars_reg.Gam;
            % end

            S = obj.S_corner{ind};
            % Phi = obj.weight_sum(Phi, parl);
            % Gam = obj.weight_sum(Gam, parl);

            sys =obj.fetch_model(S, Phi, Gam);
            % sys = cell(obj.Nss, 1);
            % for i = 1:obj.Nss
            %     sys{i} =obj.fetch_model(obj.S{i}, Phi{i}, Gam{i});
            % end
            

        end



        function plant_model = connect_model(obj, plant, ind, rho)
            %connect the model (nominal regulator equation)

            if nargin < 3
                rho = 1;
            end

            % if ~iscell(plant)
            %     plant0 = plant;
            %     plant = cell(obj.Nss, 1);
            %     for i = 1:obj.Nss
            %         plant{i} = plant0;
            %     end
            % end


            %form the model interconnection
            I = eye(obj.Nss);
            % for i = 1:obj.Nss
                parl = I(ind, :);
                model = obj.get_model(parl);
                model_rho = rhotrafo(model, rho);
                plant_model = lft(plant, model_rho);
            % end

        end
        
    end
end

