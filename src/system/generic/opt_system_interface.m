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
        
        
        discount = true; %is the subsystem exponentially discounted?
        %for a 3-mode system: could be [true, false, false]
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

        %% formation of the plant
        function [alg_psi, iqc_op, alg_loop] = build_plant_single(obj, alg, iqc_data, rho)
            %BUILD_PLANT_SINGLE build a single plant (in a switched
            %system) based on filtering the exponentially-discounted plant 
            %by an IQC
            %
            %Input:
            %   alg:        original algorithm or network
            %   iqc_data:   IQCs for the oracle uncertainties
            %   rho:        exponential discount factor
            %
            %Output:
            %   alg_psi:    filtered algorithm 
            %   iqc_op:     IQCS for the oracle uncertainties (altogether)
            %   alg_loop:   the discounted algorithm before applying the 
            %               dynamical filter (for debugging)

            %sort oracles based on the bind (exposure of repeated
            %nonlinearities)
            nop = length(obj.bind);
            [~, perm] = sort(obj.bind);

            %TODO: different number of channels per oracle
            c = obj.op{1}.c;
            P = eye(nop);
            P(:, perm) = P;
            P = kron(P, eye(c));

            wshift = c*nop;

            w_offset = ssize(alg.B, 2) - wshift;
            z_offset = ssize(alg.C, 1) - wshift;
            


            Pwp = blkdiag(P', eye(w_offset));
            Pzp = blkdiag(P, eye(z_offset));

            alg_perm = Pzp * alg * Pwp;

            %identify and get rid of the same (m=L) oracles   
            %use an explicit substitution w = m z rather than w \in F(z)
            ind_same = iqc_data.ind_same;
            
            ind_diff = setdiff(1:obj.P.nw, ind_same);
            Pd = eye(nop*c);
            Pd(:, [ind_same, ind_diff]) = Pd;
            n_same = length(ind_same);

            Pwp2 = blkdiag(Pd', eye(w_offset));
            Pzp2 = blkdiag(Pd, eye(z_offset));
            alg_perm_same = Pzp2 * alg_perm * Pwp2;

            alg_perm_m = lft(iqc_data.m_same, alg_perm_same, n_same, n_same);

            %exponentially weight the algorithm by the rate rho
            alg_rho = alg_perm_m;            
            alg_rho.A = (rho^(-1)) * alg_perm_m.A;
            alg_rho.B = (rho^(-1)) * alg_perm_m.B;
            
            %now apply the IQC to the exponentially-weighted system

            %get the iqcs for the operators
            %no loop transformations in performance
            iqc_op = iqc_data.iqc;

            loop = iqc_op.loop;
            nloop = length(loop)/2;
            alg_loop = lft(loop, alg_rho, nloop, nloop);




            %TODO: division between analysis and synthesis

            if strcmp(iqc_data.task, 'analysis')
                %form the system
                I = ss(eye(size(alg_loop.D, 2)));
                GI = [alg_loop; I];


                Psi1 = iqc_op.Psi1;
                Psi2 = iqc_op.Psi2;
                I_zp = eye(obj.P.nzp);
                I_wp = eye(obj.P.nwp);
    
                psi = blkdiag(Psi1, I_zp, Psi2, I_wp);
                
    
                alg_psi = psi * GI; 
            else
                %TODO: implement synthesis                

                n = obj.P.dump_dim;
                n.nw = n.nw - length(ind_same);
                n.nz = n.nz - length(ind_same);
                alg_psi = iqc_op.wrap_synth(alg_loop, n);
            end            
        end

        %% Dimension Counters
        function nss = Nss(obj)
            %NSS: number of subsystems
            nss = 1;
        end

        function dimn = nu(obj)
            %nu: number of states in network
            dimn = obj.P.nu;
        end

        function dimn = ny(obj)
            %nu: number of states in network
            dimn = obj.P.ny;
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
        
        function pow = discount_schedule(ordermax)
            %DISCOUNT_SCHEDULE exponential weights encountered when
            %applying the FIR filters
            %
            %[0; 1 ; 2] -> rho.^[0; 1; 2] for uniform exponential stability

            %This becomes relevant when performing shuffled systems
            %(override on switched systems) 
            %
            %TODO: switched systems
            

            if obj.discount
                pow = -(0:ordermax);
            end


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


        build_plant()
    
    end
end

