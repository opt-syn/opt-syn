classdef op_gen  < operator_interface
    %OP_GEN a general operator, that may not be the subdifferential of a
    % function in SmL    
    
    
    properties
        LINEAR = false;             %is the operator a linear map?
        prop = {'monotone', 0};     %properties of the operator
                                    %monotone
                                    %cocoercive
                                    %lipschitz
                                    %inverse lipschitz              
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
            %PROP_COUNT: count the number of properties
              pc = size(obj.prop, 1);
        end
             
        function mu = get_mu(obj)
            %GET_MU get the (strong) monotonicity parameter
            mu = [];
            for i = 1:obj.prop_count
                if strcmp(obj.prop{i, 1}, 'monotone')
                    mu = obj.prop{i, 2};
                end
            end
        end

        function loop_out = build_loop(obj, reps)
            %BUILD_LOOP construct the loop transformation  

            mu = obj.get_mu();

            if isempty(mu)
                mu = 0;
            end
        
            loop_out = [zeros(reps), eye(reps); eye(reps), mu*eye(reps)];
            
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

        function cost = build_cost(obj, var_curr)
            %BUILD_COST create the matrices M and X


            %IQC: [x; g] with g \in F(x)
            %
            %
            sz = ssize(var_curr{1}, 1);                                 

            zz = zeros(sz);
            pc = obj.prop_count;

            mu = obj.get_mu();

            %loop transformation for the strong monotonicity (?)
            if isempty(mu)
                loop_mat = eye(sz*2);
            else
                loop_mat = [eye(sz), zz; mu*eye(sz), eye(sz)];
                % loop_mat = inv([eye(sz), zz; -mu * eye(sz), eye(sz)]);
            end

            cost = zeros(2*sz);
            
            
            for p = 1:pc
                cDBM = var_curr{p};
                csym = (cDBM+cDBM');
                    switch obj.prop{p, 1}
                        case 'monotone'
                            % mu = obj.prop{p, 2};
                            % M12 = M12 + cDBM;
                            Mcurr = [zz, cDBM; cDBM', zz];
                            Mloop = Mcurr;
                            % M22 = M22 - (cDBM + cDBM') * (mu); 
                            % 
                        case 'cocoercive'
                            beta = obj.prop{p, 2};

                            Mcurr = [zz, cDBM; cDBM', -csym*beta];
                            Mloop = loop_mat'*Mcurr*loop_mat;
    
                        case 'lipschitz'
                            L2 = obj.prop{p, 2}^2;
                            % M11 = M11 - (cDBM + cDBM'); 
                            % M22 = M22 + cDBM * L2; 

                            Mcurr = [L2*csym, zz; zz', -csym];
                            Mloop = loop_mat'*Mcurr*loop_mat;
    
                        case 'inv_lipschitz'
                            L2 = obj.prop{p, 2}^2;
                            Mcurr = [-csym, zz; zz', L2*csym];
                                                        
                            Mloop = loop_mat'*Mcurr*loop_mat;
                        
                        otherwise
                            error('op_gen: unsupported property of operator')
                    end
                    cost = cost + Mloop;
                
            end

            % cost = [M11, M12; M12', M22];
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

        function [iqc] = create_iqc_identity(obj, reps)
            %CREATE_VARS form the IQC for the general operator
            %identity IQC in psi

            %Input:             
            %   rep:    number of repetitions of the operator (non-frugal)
            %
            %Output:
            %   vars:   variables of the problem
            %   cons:   constraints in the problem (in terms of the
            %           variables directly)

            if nargin < 2
                reps = 1;
            end

            vars = struct;
            pc = obj.prop_count;
            vars.cX = cell(pc, 1);
            vars.cM = cell(pc, 1);
            
            for i = 1:pc
                %identity matrix (w.r.t zero as the other point)
                vars.cM{i} = eye(reps)/(reps*pc);
                vars.cX{i} = [];
            end
            
            

            M = obj.build_cost(vars.cM);

            X = [];

            loop = obj.build_loop(reps);            
            
            [psi1, psi2] = build_psi(obj, vars, 0, reps);

            iqc_orig = iqc_loop_split(psi1, M, loop, psi2, X);

            iqc = iqc_orig.lift(obj.c);                
            
        end


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


         
            pc = obj.prop_count;

            cM = cell(pc, 1);
            for i = 1:pc
                cM_curr = lmim(['cM_', obj.sid, '_', num2str(i)], nM, nM, 'sym');
                cM{i} = cM_curr;
            end

            % cM = lmim(['cM_', obj.sid], NM, pc, 'full');
               
            cX = cell(pc, 1);
            for i = 1:pc
                if order > 0
                    cX_curr = lmim(['cX_', obj.sid, '_', num2str(i)], nX, nX, 'sym');
                    cX{i} = cX_curr;
                end
            end
            
            vars = struct;
            vars.cM = cM;
            vars.cX = cX;
        end
    end
end

