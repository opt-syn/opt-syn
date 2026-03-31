classdef op_sim
    %OP_SIM an operator used for the purposes of simulation (algorithm
    %execution). 

    %TODO: maybe this should be placed in the operator proper?
    
    properties
        
        fw= @(z, param) []; %forward evaluation (e.g. gradient)
        bw = @(z, param) []; %backward evaluation (e.g. proximal operator)        
        f = @(z, param) [];  %function value (or function values in a game)
    end
    
    methods
        function obj = op_sim(fw, bw, f)
            %OP_SIM Construct an instance of this class
            %   operations used in the evaluation of the operator
            
            
            obj.fw = fw;
            obj.bw = bw;
            if nargin == 3
                obj.f = f;
            end
        end               
    end
end

