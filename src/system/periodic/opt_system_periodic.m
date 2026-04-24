classdef  opt_system_periodic < opt_system_switch
    %OPT_SYSTEM_PERIODIC interconnection of network and operators
    %
    %a periodic system: repeated and predictable cycle evaluation   

    
    methods
        function obj = opt_system_periodic(op, P, K, bind, tracking)
            %OPT_SYSTEM_PERIODIC constructor            
            if nargin < 4
                s = length(op);
                bind = 1:s;            
            end
            
            Nss = max(length(P), length(K));
            adj = circshift(eye(Nss), -1);
            obj@opt_system_switch(op, P, K, adj, bind)
            
        end        


        function mode_next = next_mode(obj, mode)
            %next mode in the switching sequence
            nss = obj.P.Nss;
            mode_next = 1+ mod(mode, nss);
        end
    end
end

