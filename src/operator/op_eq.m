classdef op_eq < op_quad
    %OP_EQ term used in equatlity constraint
    % Qx = b, with singular values of Q between
    %  sqrt(m) and sqrt(L).
    %
    %enforced by Q'Q x = Q' b (for Q full row rank)
    %
    % allows for noncausal multipliers
    %
   
    methods
        function obj = op_eq(smin, smax, id)
            %OP_EQ Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                id = 0;
            end
            obj@op_quad(smin^2 , smax^2, id)            

            obj.EQUALITY = true;
        end
        
    end
end

