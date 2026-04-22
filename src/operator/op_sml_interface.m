classdef op_sml_interface < operator_interface
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
        function obj = op_sml_interface(m, L, c)
            %OP_SML Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                c = 1;
            end
            obj@operator_interface(c)            

            obj.m = m;
            obj.L = L;
        end
       

        % common routines

        function se = same(obj)
            %SAME no loop transformation or IQC required
            %perfectly known oracle
            se = obj.L == obj.m;
        end
        
        function sig = sigma(obj)            
            %SIGMA used to define all IQCs
            sig = 1/(obj.L - obj.m);
        end
               

        function loop = build_loop(obj, reps)
            %BUILD_LOOP construct the loop transformation
            loop_base = [-obj.sigma, 1; 1, obj.m];

            loop = kron(loop_base, eye(reps));
        end
        
    end
end

