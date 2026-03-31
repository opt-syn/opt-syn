classdef  opt_system
    %OPT_SYSTEM interconnection of network and operators
    %TODO: may be abstracted into an interface
    
    properties
        op; %a cell of operators (op_sim for simulation, op_? for analysis/synthesis)
        P;  %network
        K;  %controller
        bind; %which operators go to which output ports        
    end
    
    methods
        function obj = opt_system(op, P, K, bind)
            %OPT_SYSTEM constructor
            obj.op = op;
            obj.P = P;
            obj.K = K;
            if nargin < 3
                s = length(obj.op);
                obj.bind = 1:s;
            else
                obj.bind = bind;
            end
        end        

        %TODO: allow for parameterized systems

        function Pcurr = get_P(obj, param)
            Pcurr = obj.P.ss();
        end

        function Kcurr = get_K(obj, param)
            Kcurr = obj.K();
        end

        function sys_alg = get_alg(obj, param)
            %close the loop of the algorithm
            sys_alg = lft(obj.get_P(param), obj.get_K(param));
        end

        function [y, u] = get_internal_signals(obj, param, x_all, w_all)
            %extract the internal signals from the interconnection (y, u) 
            %using the well-posedness expression

            %Input:
            %   x_all:      all states of network and controller
            %   w_all:      all inputs to the network (except u)
            %   y_all:      all outputs to the network

            Kcurr = obj.get_K(param);
            Pcurr = obj.get_P(param);

            [nu, ny] = size(Kcurr.D);

            DK = Kcurr.D;
            DP = Pcurr.D((end-ny):end, (end-nu):end);

            CyP = Pcurr.C((end-ny):end, :);
            D1P = Pcurr.D((end-ny):end, 1:(end-nu));
            well_posed_mat = [eye(nu), -DK;
                              -DP, eye(ny)];


            Cblock = [CyP; Kcurr.C] * x_all;
            Dblock = [D1P * w_all; Kcurr.D * 0];

            revert = well_posed_mat \ (Cblock + Dblock);


            y = revert(1:ny);
            u = revert((ny+1):end);
        end

        function dimn = n(obj)
            %n: number of states
            dimn = length(obj.K.A) + obj.P.nx;
        end

        function op_out = get_op(obj, i)
            %get the operator at index i
            op_out = obj.op{obj.bind(i)};
        end

    end
end

