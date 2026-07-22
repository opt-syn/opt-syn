classdef quad_param
    %QUAD_PARAM describe the quadratic performance specification
    %
    
    
    properties
        Q; %         
        S; %cross term
        U; %
        T; %
    end
    
    methods
        function obj = quad_param(inputArg1,inputArg2)
            %QUAD_PERF Construct an instance of this class
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

