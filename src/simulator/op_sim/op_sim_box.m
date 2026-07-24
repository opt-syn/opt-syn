classdef op_sim_box < op_sim_interface
    %OP_SIM_BOX a projection onto a box    
    %this includes a hard :math:`L_\infty` norm as a special case

    
    properties        
        BOX = 1; %size of the box       
    end
    
    methods
        function obj = op_sim_box(BOX)
            %OP_SIM_BOX Constructor for box constraint
            
            
            obj@op_sim_interface();

            if nargin 
                obj.BOX = BOX;
            end
        end          

         function w = fw(obj, k, z, param)
            %forward evaluation of the procedure oracle w = E'(E z- b) 
            %
            %Args: 
            %   k (int): time index
            %   z:       input to oracle
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   w:       the :math:`w` such that :math:`w = F(z)`


            %nondifferentiable, don't include it.
            w =  zeros(size(z));
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


            %coordinate dimensions
            c = size(v, 1)/size(D, 1);

            %skew the box in the projection
            bskew =  (sqrt(D) \ (ones(c, 1) * obj.BOX));


            z = kron(eye(dl), sqrt(D)) * ...
                clip( kron(eye(dl), sqrt(D)) \ v, -bskew, bskew);
                        
        end


        function f_out = f(obj, k, z, param)
            %primal residual for the equality constraint
            %
            %Args: 
            %   k (int): time index
            %   z:       input to oracle           
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   f_out:   f_out = norm(Ez - b)

            %do not output
            f_out = [];
        end
    end
end

