classdef opt_config
    %OPT_CONFIG configuration options for optimization algorithm analysis
    %and synthesis
    
    properties
        %numerical tolerances 

        %generic options (for both analysis and synthesis)
        
        gen = struct('solver', 'lmilab', ...
            'verbose', 0)       

        %analysis only options        
        ana = struct('normalize_margin', 0.05);

        %synthesis only options
        syn = struct('reduced_order', false, ...
            'D_mask', []);

        %numerical tolerances for solving LMIs
        tol = struct('dia', 1e-11,  ... %tolerance for acceptable solution 
            'M', 1e-7, ...              %tolerance for dissipation constraints
            'X', 1e-7, ...              %tolerance for sign/terminal cost constraints 
            'G_max', 100, ...           %upper bound on norm of storage matrix (analysis)
            'GX_max', 100, ...          %upper bound on norm of primal storage matrix (synthesis)
            'GY_max', 100,...           %upper bound on norm of dual storage matrix (synthesis)
            'input_diss', 1e-4,...           %tolerance for strict input dissipation
            'K_max', 100);              %upper bound on norm of controller state space matrices


        %bisection routines
        %bisection goes in bisect_opts
    end
    
    methods
        function obj = opt_config()
            %OPT_CONFIG Construct an instance of this class
            %   Detailed explanation goes here            
        end
        
        function verdict = LMILAB(obj)
            %is LMILAB used?
            verdict = strcmp(obj.gen.solver, 'lmilab');
        end
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

