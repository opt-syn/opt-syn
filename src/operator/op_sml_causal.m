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

        function [vars, cons]= create_vars(obj, cons, order, reps)
            %CREATE_VARS form the variables in an IQC
            %
            %Input: 
            %   order:  order of the IQC [causal, noncausal]
            %   rep:    number of repetitions of the operator
            
            if nargin < 2
                cons = [];
            end

            if nargin < 3
                order = 0;
            end
            
            if nargin < 4
                reps = 1;
            end

            %nonnegative weights for the multipliers
            c = lmim('c', (order+1)*reps, reps);

            vars = struct('c', c);

            cons = elem_nonneg(c, cons);
            
        end       
        
        function [iqc_out, vars, cons] = create_iqc(obj, cons, order, reps)
            %CREATE_IQC form the variables in an IQC            
            if nargin < 2
                cons = [];
            end

            if nargin < 3
                order = 0;
            end
            
            if nargin < 4
                reps = 1;
            end

            [vars, cons] = obj.create_vars(cons, order, reps);

            
            %build the cost
            M = kron([0, 1; 1, 0], eye(reps));
            X = 0;

            %build the filter
            Psi1 = obj.build_psi(vars, order, reps);
            Psi2 = eye(reps);
            loop = obj.get_loop();

            iqc_out = iqc_loop_split(Psi1, M, loop, Psi2, X);
        end
    end
end

