classdef op_eq_causal < op_quad_causal
    %OP_EQ_CAUSAL term used in equatlity constraint
    % Qx = b, with singular values of Q between
    %  sqrt(m) and sqrt(L).
    %
    %enforced by Q'Q x = Q' b (for Q full row rank)
    
 
    
    methods
        function obj = op_eq_causal(smin, smax, id)
            %OP_QUAD Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                id = 0;
            end
            obj@op_quad_causal(smin^2 , smax^2, id)            

            obj.EQUALITY = true;
        end
        

    end
end

