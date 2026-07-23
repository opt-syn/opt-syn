classdef bridge_pass_through < genplant
    %BRIDGE_PASS_THROUGH no network dynamics, no performance channels
    %the simplest possible direct connection from the controller to the
    %operator. 

    
    methods
        function obj = bridge_pass_through(s, c)
            %BRIDGE_PASS_THROUGH Constructor
            %
            %Args:
            %   s (int): number of operators (including repetitions in bind)
            %   c (int): kronecker coordinate lift
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
    
            obj@genplant(G0, n);

            
            [obj.P.InputName, obj.P.OutputName] = obj.P_names(s, c);
        end


        function [InputName, OutputName] = P_names(obj, s, c)
            %P_NAMES label the input and output signals
            %
            %Args:
            %   s (int): number of operators (including repetitions in bind)
            %   c (int): kronecker coordinate lift
            %
            %Returns:
            %   InputName (cell of char):  names of inputs 
            %   OutputName (cell of char): names of outputs
            InputName = cell(2*s, 1);
            OutputName = cell(2*s, 1);
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
       
    end
end

