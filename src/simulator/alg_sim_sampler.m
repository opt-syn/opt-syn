classdef alg_sim_sampler
    %alg_sim_SAMPLER sampler for algorithm simulation 
    %random generation
    
    properties
        wp = @(k, param_in) []; %performance input
        param = @(k, param_in) []; %parameters at each time (transition rule)
        x0 = @() []; %initial state
        param0 = @() []; %initial parameters              
    end
    
    methods
        function obj = alg_sim_sampler()
            %alg_sim_samplerconstruct a sampler           
        end               
        

    end
end

