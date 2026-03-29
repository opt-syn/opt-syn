classdef op_sml < operator_interface
    %OP_SML An operator which is the subdifferential of a function in SmL:
    %
    %F = partial f, where f(x) - m norm(x, 2)^2 and L norm(x, 2)^2 - f(x)
    %are both proper, convex, and closed with -Inf < m <= L < inf
    %
    %
    % TODO: generalize to matrices m and L?
    
    properties
        m;
        L;
    end
    
    methods
        function obj = op_sml(m, L, id)
            %OP_SML Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                id = 0;
            end
            obj@operator_interface(id)            

            obj.m = m;
            obj.L = L;
        end
        
        function vars_out = create_vars(obj, order, reps)
            %CREATE_VARS form the variables in an IQC
            %
            %Input: 
            %   order:  order of the IQC [causal, noncausal]
            %   rep:    number of repetitions of the operator
            if length(order)== 1
                order = [0; order];
            end
            if nargin < 3
                reps = 1;
            end
            vars_out = [];
        end

        function con = filter_constraints(obj, vars, rho)
            %TODO: include Zames-Falb constraint
            con = [];
        end
        
        function iqc_out = create_iqc(obj, order, reps)
            %CREATE_IQC form the variables in an IQC            
            if length(order)==1
                order = [0; order];
            end
            if nargin < 3
                reps = 1;
            end
            iqc_out = [];
        end
    end
end

