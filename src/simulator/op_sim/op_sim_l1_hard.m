classdef op_sim_l1_hard < op_sim_interface
    %OP_SIM_l1_hard a projection onto an L1 ball
        

    properties
        tau=1; %radius of L1 norm ball
    end
    
    methods
        function obj = op_sim_l1_hard(tau)
            %OP_SIM_L1_hard constructor for l1 norm ball
            %   operations used in the evaluation of the operator
            %
            %Args:
            %   tau: radius

          
            
            obj@op_sim_interface();

            if nargin 
                obj.tau = tau;
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
            %   w:       the w such that w = F(z)


            %nondifferentiable, don't include it.
            w =  zeros(size(z));
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


            %coordinate dimensions
            dl = obj.blocksize(v);
            

            Dskron = kron(eye(dl), sqrt(D));
            %skew the box in the projection
            z_weighted = weighted_l1_proj(Dskron \ v, obj.tau, diag(Dskron));
            z = Dskron * z_weighted; 


            %for testing
            % z2= ProjectL1Ball( v, obj.tau, 1e-8);
            % z = z2;
            % err = z - z2;
                                  
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

