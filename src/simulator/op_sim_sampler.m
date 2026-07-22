classdef op_sim_sampler
    %OP_SIM_SAMPLER sampler for algorithm simulation 
    %random generation
    
    properties
        wp = @(param_in) []; %performance input
        x0 = @() []; %initial state
        param = @(param_in) []; %parameters at each time
        param0 = @(param_in) []; %initial parameters
    end
    
    methods
        function obj = op_sim_sampler()
            %OP_SIM_SAMPLER construct a sampler
            % obj.Property1 = inputArg1 + inputArg2;
        end
        
        function obj = set.x0(obj, x0_in)
            %Set x0
            % Args:
            %   x0_in: initial condition
            obj.x0 = x0_in;
        end

        function x0_out = get.x0(obj)
            %get the initial condition x0
            if isnumeric(obj.x0)
                x0_out  = obj.x0;
            else
                x0_out =obj.x0();
            end
        end

        function obj = set.wp(obj, wp_in)
            %Set wp
            % Args:
            %   x0_in: initial condition
            obj.wp = x0_in;
        end

        function wp_out = get.wp(obj, par_in)
            %get the initial condition x0
            % Args:
            %   par_in: input parameters
            if nargin < 2
                par_in = [];
            end
            if isnumeric(obj.wp)
                wp_out  = obj.wp;
            else
                wp_out =obj.wp(par_in);
            end
        end

        

    end
end

