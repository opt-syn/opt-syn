classdef bridge_random < bridge_pass_through
    %BRIDGE_RANDOM randomly generated state-space system 
    

    
    methods
        function obj = bridge_random(nconn, s, eigset)
            %BRIDGE_RANDOM Constructor
            %
            %
            %Args:
            %   s (int): number of operators (including repetitions in bind)                      
            %   eigset (float): maximum eigenvalue of the A matrix (optional)

         

            obj@bridge_pass_through(s);

            obj.P = drss(nconn, 2*s, 2*s);
            obj.P.Ts = 1;   
            
            if nargin == 3 && nconn>0
                obj.P.A = eigset* obj.P.A / max(abs(eig(Gd.A)));
            end

            for i = 1:length(obj.P.A)
                obj.P.StateName{i} = sprintf('x%d', i);
            end

            [obj.P.InputName, obj.P.OutputName] = obj.P_names(s);
        end
        
     
    end
end

