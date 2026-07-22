classdef opt_config_gen
    %OPT_CONFIG_GEN Configuration options for analysis and synthesis
    
    properties
        solver = 'lmilab'; %which solver to use
        verbose = false; %level of output
        impose_X = true; %apply PD constraint to the storage matrix
    end
    
    methods
        function obj = opt_config_gen()
            %OPT_CONFIG_GEN constructor
        end        
    end
end

