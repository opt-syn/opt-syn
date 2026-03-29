classdef (Abstract) opt_manager_interface
    %OPT_MANAGER_INTERFACE interface for the analysis and synthesis of
    %optimization/fixed point algorithms
    
    properties
        system;
    end
    
    methods
        function obj = opt_manager_interface(op, network, bind)
            %OPT_MANAGER_INTERFACE
            obj.op = op;
            obj.network = network;
            if nargin < 3
                s = length(obj.op);
                obj.bind = 1:s;
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

