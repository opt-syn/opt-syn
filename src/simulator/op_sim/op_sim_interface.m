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
        
    end
    
    methods(Abstract)
        fw; %forward evaluation
        bw; %backward evaluation
        f;  %function evaluation
    end
end

