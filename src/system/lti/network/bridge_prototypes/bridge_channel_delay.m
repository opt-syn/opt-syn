classdef bridge_channel_delay < bridge_pass_through
    %BRIDGE_DELAYS channel delays before and after the oracle
    

    
    methods
        function obj = bridge_channel_delay(d1, d2, c)
            %BRIDGE_CHANNEL_DELAY Construct an instance of this class
            %   Detailed explanation goes here

            if nargin < 3
                c = 1;
            end

            s = length(d1);
            obj@bridge_pass_through(s, c);

            
            G1 = ss(zeros(s*c)); %from u to z
            G2 = ss(zeros(s*c)); %from w to y   
            
            z = tf('z', 1);
            for i = 1:s
                for j = 1:c

                    ic = (i-1)*c;
                    G1(ic + j, ic + j) =  z^(-d1(i));
                    G2(ic + j, ic + j) =  z^(-d2(i));
                end
            end
        
            %get the final model
            obj.P =  obj.P * blkdiag(G2, G1);
        
            for i = 1:length(obj.P.A)
                obj.P.StateName{i} = sprintf('x%d', i);
            end

            [obj.P.InputName, obj.P.OutputName] = obj.P_names(s, c);
        end


        function [InputName, OutputName] = P_names(obj, s, c)

            % if nargin 
            InputName = cell(s*c, 1);
            OutputName = cell(s*c, 1);
            for i = 1:(2*s)
                for j = 1:c
                    ic = (i-1)*c;
                    if i <= s
                        InputName{ic+j} = sprintf('w%d_%d', i, j);
                        OutputName{ic+j} = sprintf('z%d_%d', i, j);
                    else
                        InputName{ic+j} = sprintf('u%d_%d', i-s, j);
                        OutputName{ic+j} = sprintf('y%d_%d', i-(s), j);
                    end
                end
            end
        end

    end
end

