classdef  opt_system_interface
    %OPT_SYSTEM interconnection of network and operators
    %by default is LTI (linear time invariant)
    
    properties
        op; %a cell of operators (op_sim for simulation, op_? for analysis/synthesis)
        P;  %network
        K;  %controller
        bind; %which operators go to which output ports            
        tracking; %tracking of optimal solution (struct (S, R) by default)
                  %tracking of varying gradients requires LPV/periodic/switched
                  % methods, is a TODO
        type=[];   %type of system: (e.g. lti, periodic, switched, mjls, lpv)
    end
    
    methods
        function obj = opt_system_interface(op, P, K, bind, tracking)
            %OPT_SYSTEM constructor for the system
            obj.op = op;
            obj.P = P;
            obj.K = K;
            if nargin < 4
                s = length(obj.op);
                obj.bind = 1:s;
            else
                obj.bind = bind;
            end

            %assign identifiers to the operators
            for i = 1:length(op)
                obj.op{i}.id = i;
            end
            
            if nargin >= 5
                obj.tracking = tracking;
            end
            
        end    

        %% Dimension Counters
        function nss = Nss(obj)
            %NSS: number of subsystems
            nss = 1;
        end


        function dimn = n(obj)
            %n: number of states
            dimn = obj.nxn() + obj.nxi();
        end

        function tp = get_type(obj)
            tp = obj.type;
        end
        

        %% getters
        
        %must define get_P, get_K

        function sys_alg = get_alg(obj, param)
            %close the loop of the algorithm
            if nargin < 2
                param = [];
            end
            Pcurr = obj.get_P(param);
            Kcurr = obj.get_K(param);
            sys_alg = lft(Pcurr, Kcurr);
        end

        function op_out = get_op(obj, i)
            %get the operator at index i
            op_out = obj.op{obj.bind(i)};
        end
        
        %% for simulation

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
            DP = Pcurr.D((end-ny+1):end, (end-nu+1):end);


            nxi = length(Kcurr.A);
            nx = length(Pcurr.A);
            CyP = Pcurr.C((end-ny+1):end, :);            
            D21P = Pcurr.D((end-ny+1):end, 1:(end-nu));
            well_posed_mat = [eye(nu), -DP;
                              -DK, eye(ny)];
            

            nx = size(Pcurr.A, 1);
            nxi = size(Kcurr.A, 1);
            xN = x_all(1:(nx), :);
            xi = x_all((end-nxi+1):end, :);

            Cx = CyP*xN;
            Dw = D21P*w_all;
            Cxi = Kcurr.C * xi;

            sig_rhs = [Cx + Dw; Cxi];
            revert = well_posed_mat \ sig_rhs;


            y = revert(1:ny, :);
            u = revert((ny+1):end, :);
        end


        function [Sbeta, Rbeta] = get_tracked_opt(obj)
            %GET_TRACKED_OPT get the tracked position of the optimal
            %solution
            if isempty(obj.tracking)
                Sbeta = 1;
                Rbeta = 1;
            else
                Sbeta = obj.tracking.Sbeta;
                Rbeta = obj.tracking.Rbeta;
            end
        end


        function mode_next = next_mode(obj, mode)
            %next mode in switching

            %TODO: is this actually used?
            mode_next = 1;
        end

 

        function N = get_consensus(obj, op, bind)
            %GET_CONSENSUS create the consensus matrix
            %for the regulation condition

            %which operators are equality constraints
            % EQ = cellfun(@(e) e.EQUALITY, op);
            nop = length(op);
            EQ = zeros(1, nop, 'logical');
            for i = 1:nop
                EQ(i) = op{i}.EQUALITY;
            end

            
            if all(~EQ)
                s = length(op);
                N0 = [eye(s-1); -ones(1, s-1)];
                
            else  
                s = sum(~EQ);
                N0 = full(sparse(1:s, find(~EQ), ones(s, 1), nop, s));
                % N0 = [eye(s)];
            end

            %index based on the bind 
            nbind = length(bind);
            Bind = full(sparse(1:nbind, bind, ~EQ(bind), nbind, nop));


            N = Bind * N0;
        end                

    end

    methods (Abstract)
        get_K(obj, param) %get the controller at the current parameter values
        get_P(obj, param) %get the network at the current parameter values

        nxn(obj)    %number of network states
        nxi(obj)    %number of controller states

        %regulator equations
        form_internal_model(obj)
        check_regulator(obj)
    end
end

