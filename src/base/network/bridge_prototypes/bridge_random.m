classdef bridge_random < bridge_pass_through
    %BRIDGE_RANDOM randomly generated state-space system 
    

    
    methods
        function obj = bridge_random(nconn, s, eigset)
            %BRIDGE_RANDOM Construct an instance of this class
            %   Detailed explanation goes here
         

            obj@bridge_pass_through(s);

            obj.P = drss(nconn, 2*s, 2*s);
            obj.P.Ts = 1;   
            
            if nargin == 3 && nconn>0
                obj.P.A = eigset* obj.P.A / max(abs(eig(Gd.A)));
            end

            for i = 1:length(obj.P.A)
                obj.P.StateName{i} = sprintf('x%d', i);
            end

            % [obj.P.InputName, obj.P.OutputName] = obj.P_names(s);
        end


        function [InputName, OutputName] = P_names(obj, s)
            InputName = cell(s, 1);
            OutputName = cell(s, 1);
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
        
     
    end
end

