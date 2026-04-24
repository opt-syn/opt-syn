classdef lmi_dispatch_lti < lmi_dispatch_interface
    %LMI_DISPATCH_LTI analysis and synthesis LMIs for the
    %linear-time-invariant (LTI) networks and controllers
    
    properties
        Property1
    end
    
    methods
        function obj = lmi_dispatch_lti(inputArg1,inputArg2)
            %LMI_DISPATCH_LTI Construct an instance of this class
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

