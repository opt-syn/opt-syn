classdef  opt_system_periodic_orbit < opt_system_periodic
    %OPT_SYSTEM_PERIODIC_orbit interconnection of network and operators
    %
    %a periodic system: repeated and predictable cycle evaluation  
    %
    %
    %a known matrix R encodes the periodicity (to be improved)

    properties
        R = 0;
    end
    
    methods
        function obj = opt_system_periodic_orbit(R, op, P, K, bind, tracking)
            %OPT_SYSTEM_PERIODIC constructor            
            if nargin < 5
                s = length(op);
                bind = 1:s;            
            end

            if nargin < 6
                 tracking = [];
            end

            obj@opt_system_periodic(op, P, K, bind, tracking)
            obj.R = R;
            obj.type = 'periodic_orbit';
        end        

    end
end

