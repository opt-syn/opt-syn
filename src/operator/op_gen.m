classdef op_gen  < operator_interface
    %OP_GEN a general operator set-valued map
    
    
    properties
        monotone = []; %monotonicity constant
        cocoercive =[]; %cocoercive constant
        lipschitz = []; %lipschitz  constant
        inverse_lipschitz=[]; %inverse_lipschitz constant
    end
    
    methods
        function obj = op_gen(c)
            %OP_GEN Construct a general operator (possibly a set-valued
            %map that does not have a potential function) 
            %
            %fill in properties by .set assignments after the constructor
            %
            %Args:            
            %   c:  dimension of coordinate lift


            if nargin < 1
                c = 1;
            end
            obj@operator_interface(c)   
                        
            obj.c = c;


        end        

        function prop = prop_report(obj)
            %PROP_REPORT get list of all properties

            prop_names = {'monotone', 'lipschitz', ...
                'inverse_lipschitz', 'cocoercive'};

            prop= {};
            for i = 1:length(prop_names)
                curr = obj.(prop_names{i});

                if ~isempty(curr) && isnumeric(curr)
                    prop_new = {prop_names{i}, curr};
                    prop = [prop; prop_new];
                end
            end


        end

        function pc = prop_count(obj)
            %PROP_COUNT count the number of properties     
            %
            %Returns:
            %   pc:     number of properties
            pc= size(obj.prop_report(), 1);
        end
             
        function mu = get_mu(obj)
            %GET_MU get the (strong) monotonicity parameter
            %
            %Returns:
            %   mu:     strong monotonicity parameter
            
            mu = obj.monotone;
            
        end

        function loop_out = build_loop(obj, reps)
            %BUILD_LOOP construct the signal transformation matrix
            %
            %Args:
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   loop_out: signal transformation matrix for the operator



            mu = obj.get_mu();


            if isempty(mu)
                beta = obj.cocoercive;
                if isempty(beta)
                    loop_out = eye(2*reps);
                else
                    loop_out = [eye(reps), -beta*eye(reps); zeros(reps), eye(reps)];
                end

            else
                loop_out = [eye(reps), zeros(reps); -mu*eye(reps), eye(reps)];
            end
            

            % loop_out = [zeros(reps), eye(reps); eye(reps), mu*eye(reps)];
            
        end

        function [psi1, psi2] = build_psi(obj, vars, order, reps)
            %BUILD_PSI construct the filter for the general operator                        
            %
            %Args:
            %   vars:   variables of the problem    
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   psi1: filter on output (causal)
            %   psi2: filter on input (noncausal components)

            if ~isscalar(order)
                order = sum(order);
            end
            
            psi_base = obj.build_psi_fir(order, reps);

            psi1 = psi_base;
            psi2 = psi_base;
        end

        function cs = csum_psi(obj, vars)
            %A proxy to normalize the filter coefficients, reducing degrees
            %of freedom in the Analysis problem
            %
            %Args:
            %   vars:   variables of the problem 
            %Returns:
            %   cs: the sum of nonnegative variables
            
            % cs = 1;

            % [n, m] = dim(vars.cM);
            % cs = ones(1, n) * vars.cM * ones(m, 1);

            cs = 0;
            for i = 1:obj.prop_count
                cs = cs + trace(vars.cM{i}) + trace(vars.cX{i});
            end
            % cs = 1;

        end

        function psi = build_psi_fir(obj, order, reps)
            %BUILD_PSI_FIR form the fir filter [1; z^-1; z^-2; z^-3..] 
            % repeated by repetitions in reps 
            %
            %Args:
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %Returns:
            %   psi: the filter in the IQC

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
            %
            %Args:
            %   vars:   variables of the problem 
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %Returns:
            %   M_out: the running cost
            
            M_out = obj.build_cost(vars.cM);
        end

        function X_out = build_X(obj, vars, order, reps)
            %BUILD_X create the terminal cost X
            %
            %Args:
            %   vars:   variables of the problem 
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %Returns:
            %   X_out: the terminal cost

            if order > 0
                % X_out = obj.build_cost(order*reps, vars.cX);
                X_out = obj.build_cost(vars.cX);
            else
                X_out = [];
            end
        end

        function cost = build_cost(obj, var_curr)
            %BUILD_COST create the matrices M and X
            %Args:
            %   vars_curr: current variables (for M or X)            
            %Returns:
            %   cost: the cost matrix (M or X)

            %IQC: [g; x] with g \in F(x)
            %
            %
            sz = ssize(var_curr{1}, 1);                                 

            zz = zeros(sz);
            
            loop_mat = inv(obj.build_loop(sz));
            % mu = obj.get_mu();
            % 
            % %loop transformation for the strong monotonicity (?)
            % %make sure that 0 is in the set of considered uncertainties
            % if isempty(mu)
            %     loop_mat = eye(sz*2);
            % else
            %     % loop_mat = kron([mu, 1; 1, 0], eye(sz));
            %     loop_mat = kron([1, 0; mu, 1], eye(sz));                
            % end

            cost = zeros(2*sz);
            
            prop = obj.prop_report();
            pc = obj.prop_count;
            for p = 1:pc
                cDBM = var_curr{p};
                csym = (cDBM+cDBM');
                    switch prop{p, 1}
                        case 'monotone'
                            % mu = prop{p, 2};
                            % M12 = M12 + cDBM;
                            Mcurr = [zz, cDBM; cDBM', zz];
                            Mloop = Mcurr;


                            % M22 = M22 - (cDBM + cDBM') * (mu); 
                            % 
                        case 'cocoercive'
                            beta = prop{p, 2};

                            Mcurr = [zz, cDBM'; cDBM, -csym*(beta)];
                            Mloop = loop_mat'*Mcurr*loop_mat;
    
                        case 'lipschitz'
                            L2 = prop{p, 2}^2;
                            % M11 = M11 - (cDBM + cDBM'); 
                            % M22 = M22 + cDBM * L2; 

                            Mcurr = [L2*csym, zz; zz', -csym];
                            Mloop = loop_mat'*Mcurr*loop_mat;
    
                        case 'inv_lipschitz'
                            L2 = prop{p, 2}^2;
                            Mcurr = [-csym, zz; zz', L2*csym];
                                                        
                            Mloop = loop_mat'*Mcurr*loop_mat;
                        
                        otherwise
                            error('op_gen: unsupported property of operator')
                    end
                    cost = cost + Mloop;
                
            end

 
        end


        function cons = filter_constraints(obj, cons, order, vars, rho_sched, iqc_out)
            %constraints on the filter coefficients (variables)
            %
            %Args:
            %   cons:   accumulated constraints
            %   vars:   variables of the problem             
            %   rho_sched:  which times should be discounted
            %   iqc_out:    the IQC under consideration            
            %Returns:
            %   cons:   accumulated constraints

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
                %impose DHD constraints
                [cons] = dhd_impose(vars.cM{i}, cons, obj.LMILAB);
                if order > 0
                    [cons] = dhd_impose(vars.cX{i}, cons, obj.LMILAB);
                end
            
            end

        end

        function [iqc] = create_iqc_identity(obj, reps)
            %CREATE_IQC_IDENTITY form a valid IQC satisfied by the general
            %operator. This is used as a warm start in synthesis.
            %
            %Args:             
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   iqc (iqc_loop_split): a valid IQC with no dynamics    

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
            %Args: 
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   vars:   variables of the problem            


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

