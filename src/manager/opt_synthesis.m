classdef opt_synthesis < opt_manager_interface
    %OPT_SYNTHESIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        reduced_order = 0;
    end
    
    methods
        function obj = opt_synthesis(sys)
            %OPT_SYNTHESIS Construct an instance of this class
            %   Detailed explanation goes here

            obj@opt_manager_interface(sys);
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

