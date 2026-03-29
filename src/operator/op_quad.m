classdef op_quad < op_sml
    %OP_QUAD Summary of this class goes here
    %   Detailed explanation goes here
    
 
    
    methods
        function obj = op_quad(m, L, id)
            %OP_QUAD Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                id = 0;
            end
            obj@op_sml(m , L, id)            

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
            %TODO: include Positive-Real constraint
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

