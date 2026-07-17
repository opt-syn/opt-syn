classdef bisect_opts < handle
    %SWEEP_OPTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        val_range = [1e-4, 2];
        bisect_rho = true;
        spec_ind = 1;
        val_init = [];
        tol = 1e-4; %tolerance for upper bound - lower bound in bisection
        warm_factor = 0.1; %if an initial guess is supplied, expand for an 
        % initial range
        Niter = 4; %number of alternating iterations

        backoff = 1e-4;
        bisect = true;

    end
    
    methods
        function obj = bisect_opts()
            %BISECT_OPTS Construct an instance of this class
            %   Detailed explanation goes here
            % obj.Property1 = inputArg1 + inputArg2;
        end
        

    end
end

