classdef alg_sim_sampler
    %alg_sim_SAMPLER sampler for algorithm simulation 
    %random generation
    
    properties
        wp = @(param_in) []; %performance input
        x0 = @() []; %initial state
        param = @(param_in) []; %parameters at each time (transition rule)
        param0 = @(param_in) []; %initial parameters
    end
    
    methods
        function obj = alg_sim_sampler()
            %alg_sim_samplerconstruct a sampler           
        end               
        

    end
end

