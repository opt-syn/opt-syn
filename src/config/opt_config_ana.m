classdef opt_config_ana
    %OPT_CONFIG_ANA Configuration options for analysis
    
    properties
        normalize_margin = 0.05; %The sum of the coefficients in the analysis program should be within the interval [s-margin, s+margin].
    end
    
    methods
        function obj = opt_config_ana()
            %OPT_CONFIG_ANA constructor
        end        
    end
end

