classdef  opt_system_periodic < opt_system_switch
    %OPT_SYSTEM_PERIODIC interconnection of network and operators
    %polytopic setting: a cell A = sum theta_i A_i for parameters theta_i
    %
    %useful for switched systems, periodic systems, and LPV systems

    
    methods
        function obj = opt_system_periodic(op, P, K, bind)
            %OPT_SYSTEM_PERIODIC constructor            
            if nargin < 4
                s = length(op);
                bind = 1:s;            
            end
            
            obj@opt_system_switch(op, P, K, G, bind)
            obj.G = circshift(eye(obj.Nss), -1);
        end        
    end
end

