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
            %OP_GEN Construct an instance of this class
            %   Detailed explanation goes here
            
        
            
             if nargin < 2
                id = 0;
            end
            obj@operator_interface(id)   
            
            if nargin > 0
                obj.prop = prop;
            end


        end

        function pc = prop_count(obj)
              pc = size(obj.prop, 1);
        end
        
     

        function [iqc, vars, cons] = create_iqc(obj, cons, order, reps)
            %CREATE_VARS form the IQC for the general operator

            %Input: 
            %   order:  order of the IQC [number of lags]
            %   rep:    number of repetitions of the operator (non-frugal)
            %
            %Output:
            %   vars:   variables of the problem
            %   cons:   constraints in the problem (in terms of the
            %           variables directly)

            [vars, cons] = obj.create_vars(cons, order, reps);


    
            %form the IQC

            M = obj.build_cost((order+1)* reps, vars.cM);
            X = obj.build_cost(order * reps, vars.cX);

            psi_base = obj.build_psi_fir(order, reps);

            psi1 = psi_base;
            psi2 = psi_base;

            orep = order*reps;
            loop = [zeros(orep), eye(orep); eye(orep), zeros(orep)];


            iqc = iqc_loop_split(psi1, M, loop, psi2, X);

        end

        function psi = build_psi_fir(obj, order, reps)
            %BUILD_PSI_FIR form the fir filter [1; z^-1; z^-2; z^-3..] 
            % repeated by reps 
            [Af, Bf] = block_fir(order);
            Cf = Af;
            Df = Bf;

            Af_all = kron(eye(reps), Af);
            Bf_all = kron(eye(reps), Bf);
            Cf_all = kron(eye(reps), Cf);
            Df_all = kron(eye(reps), Df);

            psi = ss(Af_all, Bf_all, Cf_all, Df_all, 1);

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

        function [vars, cons] = create_vars(obj, cons, order, reps)
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
                cons = [];
            end

            if nargin < 3
                order = 0;
            end
            if nargin < 4
                reps = 1;
            end

     
            %declare the variables
            n = (order+1) * reps;
            N = n + n*(n-1)/2;


            nX = order * reps;
            NX = nX + nX*(nX-1)/2;

            pc = obj.prop_count;


            cM = lmim('cM', N, pc, 'full');
            cX= lmim('cX', NX, pc, 'full');
            
            vars = struct('cM', cM, 'cX', cX);

            %declare the constraints            
            cons = elem_nonneg(cM, cons);
            cons = elem_nonneg(cX, cons);          

        end
    end
end

