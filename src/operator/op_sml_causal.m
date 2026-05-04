classdef op_sml_causal_lam < op_sml_causal
    %OP_SML_CAUSAL_lam An operator which is the subdifferential of a function in SmL:
    %
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
        function obj = op_sml_causal_lam(m, L, c)
            %OP_SML Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                c = 1;
            end
            obj@op_sml_causal(m ,L, c)            
        end



        function M = build_M(obj, vars, order, reps);
            %BUILD_M create the running cost M
            M = kron([0, 1; 1, 0], eye(reps));
        end

        function X_out = build_X(obj, vars, order, reps)
            %BUILD_X create the terminal cost X
            X_out = 0;
        end

        function [iqc, vars, cons] = create_iqc(obj, cons, order, reps)
            %create the IQC
            if length(order)>1
                order = sum(order);
            end

            [iqc, vars, cons] = create_iqc@op_sml_interface(obj, cons, order, reps);
        end


        function sm = same(obj)
            %SAME: is there any uncertainty in this oracle?
            sm = (obj.m == obj.L);
        end


        function sm = get_same(obj, reps)
            %GET_SAME: is there any uncertainty in this oracle?
            sm = kron(obj.m, eye(reps));
        end

        function cs = csum_psi(obj, vars)
            cs = trace(vars.Df);
        end

        %% subsidiary creation routines
        function [vars] = create_vars(obj, order, reps)
            %CREATE_VARS form the variables in an IQC
            %
            %Input: 
            %   order:  order of the IQC [causal, noncausal]
            %   rep:    number of repetitions of the operator
 
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
            %constraints on the filter coefficients


            cons = elem_nonneg(-vars.Cf, cons, obj.LMILAB);
            r = ssize(vars.Df, 2);
            nsched = size(rho_sched, 2);
            for i = 1:nsched
                curr_sched = rho_sched(:, i);
                M1 = [vars.Cf, vars.Df]*curr_sched;
                cons = elem_nonneg(M1, cons, obj.LMILAB);

            end
            
            cons = elem_nonneg_offdiag(-vars.Df, cons, obj.LMILAB);            

        end

        function [Psi1, Psi2] = build_psi(obj, vars, order, reps)
            %BUILD_PSI construct the filter for the SML function
            %
            %use Zames-Falb multipliers to do this

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

