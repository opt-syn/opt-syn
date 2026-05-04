classdef op_gen  < operator_interface
    %OP_GEN a general operator, that may not be the subdifferential of a
    % function in SmL    
    
    
    properties
        LINEAR = false;
        prop = {'monotone', 0};
        %other properties of an operator?
        tol_cX = 1000; %tolerance for the terminal cost
    end
    
    methods
        function obj = op_gen(prop, c)
            %OP_GEN Construct a general operator (possibly a set-valued
            %map that does not have a potential function)            
            if nargin < 2
                c = 1;
            end
            obj@operator_interface(c)   
            
            if nargin > 0
                obj.prop = prop;
            end
            obj.c = c;


        end

        function pc = prop_count(obj)
              pc = size(obj.prop, 1);
        end
             
        function loop_out = build_loop(obj, reps)
            %BUILD_LOOP construct the loop transformation
            loop_out = [zeros(reps), eye(reps); eye(reps), zeros(reps)];
        end

        function [psi1, psi2] = build_psi(obj, vars, order, reps)
            %BUILD_PSI construct the filter for the general operator
            
            psi_base = obj.build_psi_fir(order, reps);

            psi1 = psi_base;
            psi2 = psi_base;
        end

        function cs = csum_psi(obj, vars)
            %A proxy to normalize the multipliers
            % cs = 1;

            % [n, m] = dim(vars.cM);
            % cs = ones(1, n) * vars.cM * ones(m, 1);

            cs = 0;
            for i = 1:obj.prop_count
                cs = cs + trace(vars.cM{i}) + trace(vars.cX{i});
            end

        end

        function psi = build_psi_fir(obj, order, reps)
            %BUILD_PSI_FIR form the fir filter [1; z^-1; z^-2; z^-3..] 
            % repeated by reps 
            if order > 0
            [Af, Bf] = block_fir(order);
            
            Cf = [Af; [zeros(1, order-1), 1]];
            Df = [Bf; 0];
            else
                Af = zeros(0, 0);
                Bf = zeros(0, 1);
                Cf = zeros(1, 0);
                Df = 1;
            end

            Af_all = kron(eye(reps), Af);
            Bf_all = kron(eye(reps), Bf);
            Cf_all = kron(eye(reps), Cf);
            Df_all = kron(eye(reps), Df);

            psi = ss(Af_all, Bf_all, Cf_all, Df_all, 1);

        end

        function M_out = build_M(obj, vars, order, reps)
            %BUILD_M create the running cost M
            M_out = obj.build_cost(vars.cM);
        end

        function X_out = build_X(obj, vars, order, reps)
            %BUILD_X create the terminal cost X
            if order > 0
                % X_out = obj.build_cost(order*reps, vars.cX);
                X_out = obj.build_cost(vars.cX);
            else
                X_out = [];
            end
        end

        function cost = build_cost(obj, sz, var_curr)
            %build the DHD basis



            sz = ssize(var_curr{1}, 1);            
                        
            M11 = zeros(sz);
            M12 = zeros(sz);
            M22 = zeros(sz);

            count = 1;
            for p = 1:pc
                cDBM = var_curr{p};
                    switch obj.prop{p, 1}
                    case 'monotone'
                        mu = obj.prop{p, 2};
                        M12 = M12 + cDBM;
                        M22 = M22 - (cDBM + cDBM') * mu/2; 
                
                    case 'cocoercive'
                        beta = obj.prop{p, 2};
                        M12 = M12 + cDBM;
                        M11 = M11 - (cDBM + cDBM') * beta/2; 

                    case 'lipschitz'
                        L2 = obj.prop{p, 2}^2;
                        M11 = M11 - (cDBM + cDBM')/2; 
                        M22 = M22 + cDBM * L2; 

                    case 'inv_lipschitz'
                        L2 = obj.prop{p, 2}^2;
                        M22 = M22 - (cDBM + cDBM')/2; 
                        M11 = M11 + cDBM * L2;                            
                end
                
            end

            cost = [M11, M12; M12', M22];
        end


        function cons = filter_constraints(obj, cons, order, vars, rho_sched, iqc_out)
            %constraints on the filter coefficients


            %TODO: with rho schedule

            % cons = elem_nonneg(vars.cM, cons);
            % cons = elem_nonneg(vars.cX, cons); 

            %discounting with M
            reps = ssize(vars.cM{1}, 1)/(order+1);           
            nsched = size(rho_sched, 2);

            %discounting with X
            rho_sched_drop = rho_sched(1:end-1, :);
            rho_sched_drop = unique(rho_sched_drop', 'rows')'; 
            nschedX = size(rho_sched_drop, 2);

            for i = 1:obj.prop_count
                %exponential discounting of M
                if obj.LINEAR
                    [cons] = dhd_impose(vars.cM{i}, cons, obj.LMILAB);
                    if order > 0
                        [cons] = dhd_impose(vars.cX{i}, cons, obj.LMILAB);
                    end
                else
                    for j =1:nsched
                        rho_1 = kron(diag(rho_sched(1:(order+1), j)), eye(reps));
                        rho_2 = rho_1;
        
                        M_rho = rho_1 * vars.cM{i} * rho_2;
                        [cons] = dhd_impose(M_rho, cons, obj.LMILAB);
                    end
    
                    %DHD imposition of X
                    %exponential discounting 
                    if order > 0                    
                        for j =1:nschedX
                            rho_1 = kron(diag(rho_sched_drop(1:(order), j)), eye(reps));
                            rho_2 = rho_1;
        
                            X_rho = rho_1 * vars.cX{i} * rho_2;
                            [cons] = dhd_impose(X_rho, cons, obj.LMILAB);
                        end
                    end
                end
            end

        end


            

            % [dx1, dx2] = dim(vars.cX);
            % cXtop = ones(1, dx1)*vars.cX*ones(dx2, 1);
            % cons = elem_nonneg(obj.tol_cX - cXtop, cons); 
        


        function [vars] = create_vars(obj, order, reps)
            %CREATE_VARS form the variables in an IQC
            %
            %Input: 
            %   order:  order of the IQC [number of lags]
            %   rep:    number of repetitions of the operator (non-frugal)
            %
            %Output:
            %   vars:   variables of the problem
            %   cons:   constraints in the problem (in terms of the
            %           variables directly)


            if nargin < 2
                order = 0;
            end

            if length(order) > 1
                order = sum(order);
            end

            if nargin < 3
                reps = 1;
            end

     
            %declare the variables
            nM = (order+1) * reps;
            nX = order * reps;
            % NM = nM + nM*(nM-1)/2;


         
            pc = obj.prop_count;

            cM = cell(pc, 1);
            for i = 1:pc
                cM_curr = lmim(['cM_', obj.sid, '_', num2str(i)], nM, nM, 'full');
                cM{i} = cM_curr;
            end

            % cM = lmim(['cM_', obj.sid], NM, pc, 'full');
               
            cX = cell(pc, 1);
            if order > 0
              cX_curr = lmim(['cX_', obj.sid, '_', num2str(i)], nX, nX, 'full');
              cX{i} = cM_curr;
            end
            
            vars = struct('cM', cM, 'cX', cX);
        end
    end
end

