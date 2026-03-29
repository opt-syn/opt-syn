classdef  opt_system
    %OPT_SYSTEM interconnection of network and operators
    %TODO: may be abstracted into an interface
    
    properties
        op; %a cell of operators (op_sim for simulation, op_? for analysis/synthesis)
        network;
        bind;
    end
    
    methods
        function obj = opt_system(op, network, bind)
            %OPT_SYSTEM constructor
            obj.op = op;
            obj.network = network;
            if nargin < 3
                s = length(obj.op);
                obj.bind = 1:s;
            else
                obj.bind = bind;
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

