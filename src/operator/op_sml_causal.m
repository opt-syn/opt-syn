classdef op_sml_causal < op_sml_interface
    %OP_SML_CAUSAL An operator which is the subdifferential of a function in SmL:
    
    %F = partial f, where f(x) - m norm(x, 2)^2 and L norm(x, 2)^2 - f(x)
    %are both proper, convex, and closed with -Inf < m <= L < inf
    %
    %
    %causal multipliers only (in accordance with the AR paper)
    %
    %
    % TODO: generalize to matrices m and L?
    %
    %
    %_lam: uses a different parameterization for the zames-falb
    %coefficients
 
    methods
        function obj = op_sml_causal(m, L, c)
            %OP_SML_Causal constructor
            if nargin < 3
                c = 1;
            end
            obj@op_sml_interface(m ,L, c)            
        end



        function M = build_M(obj, vars, order, reps);
            %BUILD_M create the running cost M
            %Args:
            %   vars:   variables of the problem 
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %Returns:
            %   M_out: the running cost
        
            M0 = [0, 1; 1, 0];

            if obj.ERGODIC && ~obj.same
                sig = 1/(obj.L - obj.m);
                Msub0 = obj.m * [1, sig; sig, sig^2] + sig*[0, 0; 0, 1];
                
                I0rep = diag([1, zeros(reps-1, 1)]);
                Msub = kron(Msub0, I0rep);

            else
                Msub = zeros(2*reps);
            end
            M = kron(M0, eye(reps)) + Msub;
        end

        function X_out = build_X(obj, vars, order, reps)
            %BUILD_X create the terminal cost X
            %
            %Args:
            %   vars:   variables of the problem 
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %Returns:
            %   X_out: the terminal cost (is 0 for causal)
            X_out = 0;
        end

        function [iqc, vars, cons] = create_iqc(obj, cons, order, reps)
            %CREATE_IQC_IDENTITY form a valid IQC satisfied by the sml
            %operator. This is used as a warm start in synthesis.
            %
            %Args:             
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   iqc (iqc_loop_split): a valid IQC with no dynamics    

            if length(order)>1
                order = sum(order);
            end

            [iqc, vars, cons] = create_iqc@op_sml_interface(obj, cons, order, reps);
        end


        function cs = csum_psi(obj, vars)
            cs = trace(vars.Df);
        end

        %% subsidiary creation routines
        function [vars] = create_vars(obj, order, reps)
            %CREATE_VARS form the variables in an IQC
            %
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
            
            if nargin < 3
                reps = 1;
            end


            %filter coefficients
            if order > 0
                Cf = lmim(['Cf_', obj.sid], reps, order*reps, 'full');
            else
                Cf = zeros(reps, 0);
            end
            Df = lmim(['Df_', obj.sid], reps, reps, 'full');

           

            vars = struct('Cf', Cf, 'Df', Df);

        end    

        function cons = filter_constraints(obj, cons, order, vars, rho_sched, iqc)
            %FILTER_CONSTRAINTS constraints on the filter coefficients            
            %Zames-Falb DHD constraints with terminal cost
            %
            %Args:
            %   cons:   accumulated constraints
            %   vars:   variables of the problem             
            %   rho_sched:  which times should be discounted
            %   iqc_out:    the IQC under consideration            
            %Returns:
            %   cons:   accumulated constraints


            %nonpositivity of non-main elements
            cons = elem_nonneg(-vars.Cf, cons, obj.LMILAB);            
            cons = elem_nonneg_offdiag(-vars.Df, cons, obj.LMILAB);            

            %positivity of main elements
            reps = ssize(vars.Df, 2);
            nsched = size(rho_sched, 2);
            for i = 1:nsched
                curr_sched = kron(rho_sched((order+1):-1:1, i), ones(reps, 1));
                M1 = [vars.Cf, vars.Df]*curr_sched;
                cons = elem_nonneg(M1, cons, obj.LMILAB);
            end

        end

        function [Psi1, Psi2] = build_psi(obj, vars, order, reps)
            %BUILD_PSI construct the zames-falb filter for the SML function
            %
            %Args:
            %   vars:   variables of the problem    
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   psi1: filter on output (causal)
            %   psi2: filter on input (noncausal components)

            [Af0, Bf0] = block_fir(order);
            Af = kron(eye(reps), Af0 );
            Bf = kron(eye(reps), Bf0);                      
            Cf = vars.Cf;
            Df = vars.Df;
            
            Psi1 = sdpss(Af, Bf, Cf, Df);
            Psi2 = ss(eye(reps));
        end
    end
end

