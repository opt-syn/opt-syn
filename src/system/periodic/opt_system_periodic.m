classdef  opt_system_periodic < opt_system_switched
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

            if nargin < 5
                 tracking = [];
            end

            Nss = max(P.Nss);
            adj = circshift(eye(Nss), -1);
            obj@opt_system_switched(op, P, K, adj, bind, tracking)
            obj.type = 'periodic';
        end        


        function [Sbeta, Rbeta] = get_tracked_opt(obj)
            %GET_TRACKED_OPT get the tracked position of the optimal
            %solution. allow for time-varying exosystems (periodic),
            %represented by a cell

            [Sbeta0, Rbeta0] = get_tracked_opt@opt_system_switched(obj);

            if iscell(Sbeta0)
                Sbeta = Sbeta0;
                Rbeta = Rbeta0;
            else
                Sbeta = cell(obj.P.Nss, 1);
                Rbeta = cell(obj.P.Nss, 1);

                for i = 1:obj.P.Nss
                    Sbeta{i} = Sbeta0;
                    Rbeta{i} = Rbeta0;
                end
            end
 
        end

        function tp = get_type(obj)
            %get the type of the switched system
            %is periodic!
            tp = obj.type;
        end

        function mode_next = next_mode(obj, mode)
            %next mode in the switching sequence
            nss = obj.P.Nss;
            mode_next = 1+ mod(mode, nss);
        end
    end
end

