classdef opt_config_bisect < handle
    %opt_config_bisect options used when applying bisection/alternation in the
    %manager class    
    
    properties
        val_range = [1e-4, 2]; %parameter range used in bisection
        bisect_rho = true; %should bisection occur on rho
        spec_ind = 1; %which specification should be bisected
        val_init = []; %initial value
        tol = 1e-4; %tolerance for upper bound - lower bound in bisection
        warm_factor = 0.1; %if an initial guess is supplied, expand for an initial range        

        backoff = 1e-4; %increase parameter by backoff when switching between analysis and synthesis. Try to ensure feasibility this way.
        bisect = true;  %use bisection or not

    end
    
    methods
        function obj = opt_config_bisect()
            %Constructor
        end
        

    end
end

