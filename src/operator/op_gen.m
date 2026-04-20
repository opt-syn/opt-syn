classdef op_gen  < operator_interface
    %OP_GEN a general operator, that may not be the subdifferential of a
    % function in SmL    
    
    
    properties
        prop = {'monotone', 0};
            % cocoercive = [];
            % monotone = [];
            % lipschitz = [];
            % inv_lipschitz = [];
        %other properties of an operator?
    end
    
    methods
        function obj = op_gen(prop, id)
            %OP_GEN Construct a general operator (possibly a set-valued
            %map that does not have a potential function)            
            
        
            
             if nargin < 2
                id = 0;
            end
            obj@operator_interface(id)   
            
            if nargin > 0
                obj.prop = prop;
            end
            obj.id = id;


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
            M_out = obj.build_cost((order+1)* reps, vars.cM);
        end

        function X_out = build_X(obj, vars, order, reps)
            %BUILD_X create the terminal cost X
            if order > 0
                X_out = obj.build_cost(order*reps, vars.cX);
            else
                X_out = [];
            end
        end

        function cost = build_cost(obj, sz, var_curr)
            %build the DHD basis


            [nM, pc] = dim(var_curr);            
            
            % DBM = dhd_basis_0(sz);
            M11 = zeros(sz);
            M12 = zeros(sz);
            M22 = zeros(sz);

            count = 1;
            for p = 1:pc
                for i = 1:sz
                    for j = 1:i
                    % DBM_curr = DBM(:, :, i);

                    %manually generate the DHD matrices
                    cMind = lmim_index(var_curr, i, p);
                    if j == i
                        E = zeros(sz, 1);
                        E(j)=1;
                        cbl = cMind;
                    else
                        E = zeros(sz, 2);
                        E(i, 1) = 1;
                        E(j, 2) = 1;
                        cbl = [1; -1]*cMind * [1, -1];
                        % cbl = [cMind, -cMind; -cMind, cMind];
                    end
                    

                    
                    
                        
                        % cMcurr = drep(cMind, sz);
                        % cMcurr = drep_scalar(cMind, sz);
                        % cDBM = smul(DBM_curr, );
                        % cDBM = cMind * DBM_curr;
    
                        cDBM = E * cbl * E';
                        switch obj.prop{p, 1}
                            case 'monotone'
                                mu = obj.prop{p, 2};
                                M12 = M12 + cDBM;
                                M22 = M22 - cDBM * mu; 
                        
                            case 'cocoercive'
                                beta_inv = 1/obj.prop{p, 2};
                                M12 = M12 + cDBM;
                                M11 = M11 - cDBM * beta_inv; 
    
                            case 'lipschitz'
                                L2 = obj.prop{p, 2}^2;
                                M11 = M11 - cDBM; 
                                M22 = M22 + cDBM * L2; 
    
                            case 'inv_lipschitz'
                                L2 = obj.prop{p, 2}^2;
                                M22 = M22 - cDBM; 
                                M11 = M11 + cDBM * L2; 
    
                        end
                    end
                end
            end

            cost = [M11, M12; M12', M22];
        end


        function cons = filter_constraints(obj, cons, order, vars, iqc_out)
            %constraints on the filter coefficients

            cons = elem_nonneg(vars.cM, cons);
            cons = elem_nonneg(vars.cX, cons); 
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
            if nargin < 3
                reps = 1;
            end

     
            %declare the variables
            nM = (order+1) * reps;
            NM = nM + nM*(nM-1)/2;


         
            pc = obj.prop_count;


            cM = lmim(['cM_', obj.sid], NM, pc, 'full');
               
          if order > 0
                nX = order * reps;
                NX = nX + nX*(nX-1)/2;

                cX= lmim(['cX_', obj.sid], NX, pc, 'full');
            else
                cX = [];
            end
            
            vars = struct('cM', cM, 'cX', cX);

            %declare the constraints            
         

        end
    end
end

