classdef op_sim_equality < op_sim_interface
    %OP_SIM_EQUALITY an operator used for the purposes of simulation (algorithm
    %execution). 

    %an affine operation enforcing an equality constraint Ez = b
    %
    %implemented as E'(Ez - b) = 0
    %with E full row rank

    
    properties        
        E= [];  %forward evaluation (e.g. gradient)        
        b = 0;  %function value (or function values in a game)        
        
    end
    
    methods
        function obj = op_sim_equality(E, b)
            %OP_SIM Construct an instance of this class
            %   operations used in the evaluation of the operator
            
            
            %create the operators
            %TODO: abstract to more general operators E (not stored in
            %memory as a matrix)

            
            %https://proximity-operator.net/proximityoperator.html
            obj@op_sim_interface();

            obj.E = E;
            if nargin < 2
                obj.b = 0;
            end
            obj.EQUALITY = true;
        end          

         function w = fw(k, z, param)
            %forward evaluation of the procedure oracle w = E'(E z- b) 
            %
            %Args: 
            %   k (int): time index
            %   z:       input to oracle
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   w:       the w such that w = F(z)


            w =  E' * (E*z - b);            
        end

        function z = bw(k, v, D, param)
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

            %weighted projection onto the affine subspace
            
            dl = size(v, 1)/size(D, 1);
            Dkron = kron(eye(dl), (D));

            h = size(obj.E, 1);
            zlam = [Dkron, obj.E'; obj.E, zeros(h)] \ [Dkron*v; obj.b];

            z = zlam(1:length(v));            
                        
        end

        function f_out = f(k, z, param)
            %primal residual for the equality constraint
            %
            %Args: 
            %   k (int): time index
            %   z:       input to oracle           
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   f_out:   f_out = norm(Ez - b)

            f_out = norm(E*z - b);
        end
    end
end

