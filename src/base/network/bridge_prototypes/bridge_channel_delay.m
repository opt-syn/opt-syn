classdef bridge_channel_delay < bridge_pass_through
    %BRIDGE_DELAYS channel delays before and after the oracle
    

    
    methods
        function obj = bridge_channel_delay(d1, d2, s)
            %BRIDGE_CHANNEL_DELAY Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 2
                s = length(d1);
            else
                if length(d1) ~= s
                    d1 = kron(d1, ones(s, 1));
                end
                if length(d2) ~= s
                    d2 = kron(d2, ones(s, 1));
                end
            end

            obj@bridge_pass_through(s);

            
            G1 = ss(zeros(s)); %from u to z
            G2 = ss(zeros(s)); %from w to y   
            
            z = tf('z', 1);
            for i = 1:s
                G1(i, i) =  z^(-d1(i));
                G2(i, i) =  z^(-d2(i));
            end
        
            %get the final model
            obj.P =  obj.P * blkdiag(G2, G1);
        
            for i = 1:length(obj.P.A)
                obj.P.StateName{i} = sprintf('x%d', i);
            end

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
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

