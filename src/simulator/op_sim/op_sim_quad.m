classdef op_sim_quad < op_sim_interface
    %OP_SIM_QUAD a quadratic function for algorithm simulation
    %
    %
    %:math:`f = (1/2) (z-z^*)' M (z-z^*)`


    
    properties
        M;     %quadratic matrix
        zstar; %critical point to the unconstrained quadratic minimization problem
    end

    methods
        function obj = op_sim_quad(M, zstar)
            %OP_SIM_QUAD Constructor form a quadratic function
            %
            % Args: 
            %   M: a symmetric matrix defining the quadratic form
            %   zstar: critical point to the unconstrained quadratic minimization problem
            %           is a function  of k.
            
            obj = obj@op_sim_interface();
            obj.M = M;

            if nargin < 2
                zstar = zeros(size(obj.M, 1), 1);
            end

            if isnumeric(zstar)
                obj.zstar = @(k) zstar;
            end
                    
    


        end

        function w = fw(obj, k, z, param)
            %forward evaluation of an oracle w = F(z) 
            %
            %Args: 
            %   k (int): time index
            %   z:       input to oracle
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   w:       the w such that w = F(z)

            w = obj.M* (z-obj.zstar(k));
        end

        function z = bw(obj, k, v, D, param)
            %backwards evaluation of an oracle, generalization of a 
            %proximal evaluation with preconditioner D                      
            %
            %Args: 
            %   k (int): time index
            %   v:       input to proximal oracle
            %   D:       prox parameter
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   z:       the z such that z = (I - D F)^(-1)(v)

            z = (obj.M + kron(eye(dl), inv(D))) \ ...
                (obj.M*obj.zstar(k) + kron(eye(size(obj.M, 1)/size(D, 1)), D) \ v);
        end

        function f_out = f(obj, k, z, param)
            %function value evaluation, if the operator has a potential
            %could also be a vector of function evaluations in a game.
            %
            %Args: 
            %   k (int): time index
            %   z:       input to oracle           
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   f_out:   f_out = f(z) if F = \partial f.

            f_out = 0.5*(z-obj.zstar(k))'*obj.M*(z-obj.zstar(k)); 
        end
        

    end
end

