classdef op_sim < op_sim_interface
    %OP_SIM an operator used for the purposes of simulation (algorithm
    %execution). 

    
    properties
        
        fw_func = @(k, z, param) [];  %forward evaluation (e.g. gradient)
        bw_func = @(k, z, D, param) []; %backward evaluation (e.g. proximal operator)        
        f_func = @(k, z, param) [];  %function value (or function values in a game)
        
    end
    
    methods
        function obj = op_sim(fw, bw, f)
            %OP_SIM operator used in algorithm simulation
            %
            %Args:
            %   fw: forward evalution
            %   bw: backward evalution
            %   f:  function evalution
            
            
            obj.fw_func = fw;
            obj.bw_func = bw;
            if nargin == 3
                obj.f_func = f;            
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


            w = obj.fw_func(k, z, param);
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
            %   z:       the z such that z = (I - D F)^(-1)(v)

            z = obj.bw_func(k, D, v, param);
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

            f_out = obj.f_func(k, z, param); 
        end
    end
end

