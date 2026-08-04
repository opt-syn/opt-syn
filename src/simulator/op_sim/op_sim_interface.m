classdef op_sim_interface
    %OP_SIM_INTERFACE functions to evaluate algorithm trajectories

    properties
        EQUALITY = false;
        c=1;
        id = 0;
    end
    
    methods
        function obj = op_sim_interface()
            %OP_SIM_INTERFACE   
            %blank constructor
        end

        function  dl = blocksize(obj, v)
            %BLOCKSIZE compute the coordinate lift/size of blocks
            %used for backward evaluations
            %
            %Args:
            %   v: point to perform backward evaluation
            %
            %Returns:
            %   dl: size of coordinate blocks
            
            d = size(v, 1);
            dl = d/obj.c;
        end
        
    end
    
    methods(Abstract)
        fw; %forward evaluation
        bw; %backward evaluation
        f;  %function evaluation
    end
end

