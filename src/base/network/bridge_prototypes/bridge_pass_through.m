classdef bridge_pass_through < bridge
    %BRIDGE_PASS_THROUGH no network dynamics, no performance channels
    

    
    methods
        function obj = bridge_pass_through(s)
            %BRIDGE_PASS_THROUGH Construct an instance of this class
            %   Detailed explanation goes here
            n = struct;

            n.nz = s;
            n.nw = s;
            n.nu = s;
            n.ny = s;

            G0 = ss([zeros(s), eye(s); eye(s), zeros(s)]);
    
            obj@bridge(G0, n);

            
            [obj.P.InputName, obj.P.OutputName] = obj.P_names(s);
        end


        function [InputName, OutputName] = P_names(obj, s)
            for i = 1:(2*s)
                if i <= s
                    InputName{i} = sprintf('w%d', i);
                    OutputName{i} = sprintf('z%d', i);
                else
                    InputName{i} = sprintf('u%d', i-s);
                    OutputName{i} = sprintf('y%d', i-(s));
                end
            end
        end
        
        function Eo = E(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Eo = obj.A;
        end
    end
end

