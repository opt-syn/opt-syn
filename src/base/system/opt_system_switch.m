classdef  opt_system_switch < opt_system
    %OPT_SYSTEM_SWITCH interconnection of network and operators
    %polytopic setting: a cell A = sum theta_i A_i for parameters theta_i
    %
    %useful for switched systems, periodic systems, and LPV systems
    
    
    properties        
        G;    %switching graph        
    end
    
    methods
        function obj = opt_system_switch(op, P, K, G, bind)
            %OPT_SYSTEM constructor            
            if nargin < 5
                s = length(op);
                bind = 1:s;            
            end

            % if ~iscell(P)
            %     P = {P};
            % end

            if ~iscell(K)
                K0 = K;
                K = cell(P.Nss, 1);
                for i = 1:P.Nss
                    K{i} = K0;
                end
            end

            obj@opt_system(op, P, K, bind)            
            obj.G = G;
        end        

        %TODO: allow for parameterized systems

        function mode_next = next_mode(obj, mode)
            if nargin < 2
                mode = 1;
            end

            %TODO: debug this
            g = obj.G(mode, :);
            gc = g/sum(g);

            gs = cumsum(gc);

            u = rand(1);

            % mode_next = [];
            mode_next = find(u <= gs, 1, 'first');
        end

        function Scurr = get_P(obj, param)
            Pcurr = obj.P{param.mode};
            Scurr = Pcurr.ss();
        end

        function Kcurr = get_K(obj, param)
            %TODO: override this with parameters
            Kcurr = obj.K{param.mode};
        end    
        

        %dimensions
        function dimn = nxn(obj)
            %nxn: number of states in network
            dimn = obj.P.nx;
        end

        function nss = Nss(obj)
            nss = obj.P.Nss;
        end
        function dimn = nxi(obj)
            %nxi: number of states in controller
            dimn = length(obj.K{1}.A);
        end

    end
end

