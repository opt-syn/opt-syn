classdef op_sim_lsq < op_sim_interface
    %OP_SIM_LSQ a least squares cost for algorithm simulation
    %
    %
    %:math:`f = (1/2) ||A z - b||_2^2`


    
    properties
        A; %matrix in least squares
        b; %vector in least squares
    end

    methods
        function obj = op_sim_lsq(A, b)
            %OP_SIM_LSQ Constructor form a quadratic function
            %
            % Args: 
            %   A: a rectangular matrix 
            %   b: the reference vector
            
            obj = obj@op_sim_interface();
            obj.A = A;

            if nargin < 2
                b = zeros(size(obj.A, 1), 1);
            end

            if isnumeric(b)
                obj.b = @(k) b;
            else
                obj.b = b;
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

            res = obj.A*z- obj.b(k);
            w = obj.A' * (res);
        end

        function z = bw(obj, k, D, v, param)
            %backwards evaluation of an oracle, generalization of a 
            %proximal evaluation with preconditioner D                      
            %
            %Args: 
            %   k (int): time index
            %   D:       prox parameter
            %   v:       input to proximal oracle
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   z:       the z such that z = (I + D F)^(-1)(v)

            dl = obj.blocksize(v);           
            Dkron = kron(eye(dl), (D));

            %basic unoptimized implementation
            ansvec = obj.A' * obj.b(k) + Dkron \ v;
            z = (obj.A' * obj.A + kron(eye(dl), inv(D))) \ ansvec;            
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

            res =  obj.A * z - obj.b(k);
            f_out = 0.5*norm(res, 2)^2;
        end
        

    end
end

