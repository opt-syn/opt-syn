classdef bridge_pass_through < bridge
    %BRIDGE_PASS_THROUGH no network dynamics, no performance channels
    

    
    methods
        function obj = bridge_pass_through(s, c)
            %BRIDGE_PASS_THROUGH Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 1
                c=1;
            end

            n = struct;

            n.s = s;
            n.nz = s*c;
            n.nw = s*c;
            n.nu = s*c;
            n.ny = s*c;

            G0 = ss([zeros(s*c), eye(s*c); eye(s*c), zeros(s*c)]);
    
            obj@bridge(G0, n);

            
            [obj.P.InputName, obj.P.OutputName] = obj.P_names(s, c);
        end


        function [InputName, OutputName] = P_names(obj, s, c)
            for i = 1:(2*s)
                for j = 1:c
                    indcurr = c*(i-1) + j;
                    if i <= s
                        InputName{indcurr } = sprintf('w%d_%d', i, j);
                        OutputName{indcurr } = sprintf('z%d_%d', i, j);
                    else
                        InputName{indcurr } = sprintf('u%d_%d', i-s, j);
                        OutputName{indcurr } = sprintf('y%d_%d', i-(s), j);
                    end
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

