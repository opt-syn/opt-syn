classdef bridge_poly < bridge
    %BRIDGE_POLY a network sitting between the oracle F and the controller K
    
    properties
        Property1
    end
    
    methods
        function obj = bridge_poly(inputArg1,inputArg2)
            %BRIDGE_POLY Construct an instance of this class
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

