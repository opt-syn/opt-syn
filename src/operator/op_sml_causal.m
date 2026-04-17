classdef op_sml_causal < operator_interface
    %OP_SML_CAUSAL An operator which is the subdifferential of a function in SmL:
    %
    %F = partial f, where f(x) - m norm(x, 2)^2 and L norm(x, 2)^2 - f(x)
    %are both proper, convex, and closed with -Inf < m <= L < inf
    %
    %
    %causal multipliers only (in accordance with the AR paper)
    %
    %
    % TODO: generalize to matrices m and L?
    
    properties
        m;
        L;
    end
    
    methods
        function obj = op_sml_causal(m, L, id)
            %OP_SML Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                id = 0;
            end
            obj@operator_interface(id)            

            obj.m = m;
            obj.L = L;
            obj.id = id;
        end

        function se = same(obj)
            %SAME no loop transformation or IQC required
            %perfectly known oracle
            se = obj.L == obj.m;
        end
        
        function sig = sigma(obj)            
            %SIGMA used to define all IQCs
            sig = 1/(obj.L - obj.m);
        end

        function loop = get_loop(obj, reps)
            %GET_LOOP loop transformation matrix
            loop_base = [-obj.sigma, 1; 1, obj.m];

            loop = kron(loop_base, eye(reps));
        end

      
        
        function [iqc_out, vars, cons] = create_iqc(obj, cons, order, reps)
            %CREATE_IQC form the valid IQC for the operator          
            if nargin < 2
                cons = [];
            end

            if nargin < 3
                order = 0;
            end
            
            if nargin < 4
                reps = 1;
            end
            
            if obj.same
                loop = obj.m * eye(reps);
                iqc_out = iqc_loop_split([], [], loop, [], []);
                vars = {};
            else
                [vars] = obj.create_vars(order, reps);
                
                %build the cost
                M = kron([0, 1; 1, 0], eye(reps));
                X = 0;
    
                %build the filter
                Psi1 = obj.build_psi(vars, order, reps);
                Psi2 = eye(reps);
                loop = obj.get_loop(reps);
    
                iqc_out = iqc_loop_split(Psi1, M, loop, Psi2, X);

                cons = obj.filter_constraints(cons, vars, iqc_out);
            end
        end

        function cons = filter_constraints(obj, cons, vars, iqc_out)
            %constraints on the filter coefficients

            cons = elem_nonneg(vars.c, cons);
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
            c = lmim(['c', obj.sid], (order+1)*reps, reps);

            vars = struct('c', c);

            
            
        end         

        function [Psi1] = build_psi(obj, vars, order, reps)
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
        end
    end
end

