classdef lmi_dispatch_interface
    %LMI_DISPATCH_INTERFACE store the LMI routines for each system type
    
    properties
        specs;
    end
    
    methods
        function obj = lmi_dispatch_interface(inputArg1,inputArg2)
            %LMI_DISPATCH_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

