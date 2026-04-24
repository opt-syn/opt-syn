classdef op_sml_causal < op_sml_interface
    %OP_SML_CAUSAL An operator which is the subdifferential of a function in SmL:
    %
    %F = partial f, where f(x) - m norm(x, 2)^2 and L norm(x, 2)^2 - f(x)
    %are both proper, convex, and closed with -Inf < m <= L < Inf
    %
    %
    %causal multipliers only (in accordance with the AR paper)
    %
    %
    % TODO: generalize to matrices m and L?
 
    %TODO: coordinate lift c

    methods
        function obj = op_sml_causal(m, L, c)
            %OP_SML Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                c =1;
            end
            obj@op_sml_interface(m ,L, c)            
        end



        function M = build_M(obj, vars, order, reps);
            %BUILD_M create the running cost M
            M = kron([0, 1; 1, 0], eye(reps));
        end

        function X_out = build_X(obj, vars, order, reps)
            %BUILD_X create the terminal cost X
            X_out = 0;
        end

        function cons = filter_constraints(obj, cons, order, vars, iqc)
            %FILTER_CONSTRAINTS constraints on the filter coefficients            

            %Zames-Falb constraints without terminal cost, polytopic
            %definition of the multiplier
            cons = elem_nonneg(vars.c, cons);
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
            if isempty(vars)
                cs = 0;
            else
                [m, d] = dim(vars.c);
                cs = ones( 1, m) * vars.c * ones(d, 1);
            end
        end

        %% subsidiary creation routines
        function [vars]= create_vars(obj, order, reps)
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

            %nonnegative weights for the multipliers
            c = lmim(['c_', obj.sid], (order+1)*reps, reps);

            vars = struct('c', c);
 
        end         

        function [Psi1, Psi2] = build_psi(obj, vars, order, reps)
            %BUILD_PSI construct the filter for the SML function
            %
            %use Zames-Falb multipliers to do so

            c = vars.c;

            [Af0, Bf0] = block_fir(order);
            Af = kron(eye(reps), Af0 );
            Bf = kron(eye(reps), Bf0);                      
            Cf = zeros(reps, order*reps);
            Df = zeros(reps);
            % Cf = [];
            % Df = [];

            %now fill in the terms

            %these are the zames-falb offsets
            %add them all up
            C_center = [zeros(reps, order*reps); -eye(order*reps)];
            D_off = [zeros(1, reps-1); -eye(reps-1); zeros(order*reps, reps-1)];
            D_main = ones((order+1)*reps, 1);

            C_right = C_center;
            D_right = [D_main, D_off];
            
            for i = 1:reps
                ccurr = lmim_index(c, [], i)';
                Cf_curr = ccurr * C_right;
                
                D_right_curr = circshift(D_right, i-1, 2);
                Df_curr = ccurr * D_right_curr;

                E_curr = zeros(reps, 1);
                E_curr(i) = 1;
                Cf = Cf + E_curr*Cf_curr;
                Df= Df + E_curr*Df_curr;
            end

            Psi1 = sdpss(Af, Bf, Cf, Df);
            Psi2 = ss(eye(reps));
        end
    end
end

