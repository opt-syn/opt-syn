classdef (Abstract) alg_sim_interface
    %ALG_SIM_INTERFACE execution of the algorithmic interconnection
    
    properties
        sys;
        blocksize;
    end
    
    methods 
        function obj = alg_sim_interface(sys, blocksize)
            %ALG_SIM_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function ssim = sim(obj, T, x0)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end

    methods (Abstract)        
        get_next_param(obj, param_old);     %
        get_next_system(obj, param_old);
        sample_input(obj, param);
    end
end

