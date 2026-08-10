classdef  opt_system_interface
    %OPT_SYSTEM interconnection of network and operators
    %by default is LTI (linear time invariant)
    
    properties
        op; %a cell of operators (op_sim for simulation, op_[] for analysis/synthesis)
        P;  %network
        K;  %controller
        bind; %which operators go to which output ports, for repeated evaluations
        type=[];   %type of system: (e.g. lti, periodic, switched, mjls, lpv)    
        discount = true; %is the subsystem exponentially discounted?        
        tracking; %tracking of optimal solution (struct (S, R) by default)
                  %tracking of varying gradients requires LPV/periodic/switched
                  % methods, is a TODO

        
        %for a 3-mode system: could be [true, false, false]
    end
    
    methods
        function obj = opt_system_interface(op, P, K, bind, tracking)
            %OPT_SYSTEM constructor for the system
            % Args:
            %   P (genplant): network
            %   K (genplant): controller
            %   bind (int array):   indices for repeated oracle evaluations
            %   tracking (struct): exosystem to track the optimal solution

            if ~iscell(op)
                op =  {op};
            end
            
            obj.op = op;


            obj.P = P;
            obj.K = K;

            %configure the bind/repetitions
            s = length(obj.op);
            if nargin < 4                
                obj.bind = 1:s;
            else
                obj.bind = bind;
            end

            %update the coordinate dimensions of the operators
            nbind = length(obj.bind);
            if isempty(P)
                if isempty(K)
                    c = 1;                    
                else
                    if iscell(K)
                        c = size(K{1}.D, 1)/nbind;
                    else
                        c = size(K.D, 1)/nbind;
                    end                    
                end
                obj.P = bridge_pass_through(nbind, c);
            else
                c = P.nz/nbind;
            end

            %iterate through operators
            for i = 1:s
                %enforce the coordinate dimension
                obj.op{i}.c = c;

                %assign identifiers to the operators
                obj.op{i}.id = i;
            end
            
            if nargin >= 5
                obj.tracking = tracking;
            end
            
        end    

        %% formation of the plant
        function [alg_psi, iqc_op, alg_loop] = build_plant_single(obj, alg, iqc_data)
            %BUILD_PLANT_SINGLE build a single plant (in a switched
            %system) based on filtering the exponentially-discounted plant 
            %by an IQC
            %
            %Args:
            %   alg:        original algorithm or network
            %   iqc_data:   IQCs for the oracle uncertainties
            %
            %Return:
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
            
            ind_diff = setdiff(1:(c*nop), ind_same);
            Pd = eye(nop*c);
            Pd(:, [ind_same, ind_diff]) = Pd;
            n_same = length(ind_same);

            Pwp2 = blkdiag(Pd', eye(w_offset));
            Pzp2 = blkdiag(Pd, eye(z_offset));
            alg_perm_same = Pzp2 * alg_perm * Pwp2;

            alg_perm_m = lft(iqc_data.m_same, alg_perm_same, n_same, n_same);

    

            %DO NOT exponentially weight the algorithm by the rate rho
            % alg_perm_m.A = (rho^(-1)) * alg_perm_m.A;
            % alg_perm_m.B = (rho^(-1)) * alg_perm_m.B;
            % 
            % instead, use rho-hard IQCs to allow for different weights in
            % different specifications
            
            %now apply the IQC to the exponentially-weighted system

            %get the iqcs for the operators
            %no loop transformations in performance
            iqc_op = iqc_data.iqc;

            if isempty(iqc_op)
                %all operators have explicit formulae
                %no uncertainty is introduced
                alg_loop = alg_perm_m;
                alg_psi = genplant(alg_perm_m);
            else
                loop = iqc_op.loop;
                nloop = length(loop)/2;
                alg_loop = lft(loop, alg_perm_m, nloop, nloop);
                
    
    
    
                %TODO: division between analysis and synthesis
    
                if strcmp(iqc_data.task, 'analysis')
                    %form the system
                    I = ss(eye(ssize(alg_loop.D, 2)));
                    GI = [alg_loop; I];
    
                    
    
                    Psi1 = iqc_op.Psi1;
                    Psi2 = iqc_op.Psi2;
                    I_zp = eye(obj.P.nzp);
                    I_wp = eye(obj.P.nwp);
        
                    psi = blkdiag(Psi1, I_zp, Psi2, I_wp);
                    
        

                    %ordering of channels
                    %[p, zp, q, wp]                    
                    Gpsi = psi * GI;

                    %then perform a reordering
                    %[p, q, wp, zp]
                    %this makes the performance calls easier
                    %
                    %because performance is stored in supply_quad (Schur), 
                    %while the M matrix is generally not in Schur-complement form

                    n1 = iqc_op.np;
                    n2 = iqc_op.nq;
                    nzp = obj.P.nzp;
                    nwp = obj.P.nwp;

                    i_order = 1:(n1+n2+nzp+nwp);
                    j_order = [1:n1, (n1+nzp)+(1:n2), (n1+nzp+n2) + (1:nwp), n1 + (1:nzp)];
                    v_order = ones(n1+n2+nzp+nwp, 1);
                    P_order = full(sparse(i_order, j_order, v_order));

                    %permute the order of the channels
                    Gpsi_perm = P_order * Gpsi;


                    %these dimension counts are very inelegant
                    alg_psi = genplant(Gpsi_perm); 
                    alg_psi.nz = n1 + n2;
                    alg_psi.nw = obj.P.nw;
                    alg_psi.nzp = obj.P.nzp + obj.P.nwp;
                    alg_psi.nwp = obj.P.nwp;
                else

                    n = obj.P.dump_dim;
                    
                    %special case for reduced-order control + performance
                    if ~iqc_data.augmented
                        %standard control
                        n.nw = n.nw - length(ind_same);
                        n.nz = n.nz - length(ind_same);

                        alg_psi = iqc_op.wrap_synth(alg_loop, n);
                    else
                        %reduced-order control: pop the performance
                        %channels next to the regulation channels
                        nwp_orig = n.nw - c*nop;
                        nzp_orig = n.nz - c*nop;

                        norig = n;

                        n.nw = n.nw - nwp_orig - length(ind_same);
                        n.nz = n.nz - nzp_orig - length(ind_same);
                        n.nwp = n.nwp + nwp_orig;
                        n.nzp = n.nzp + nzp_orig;


                        alg_psi = iqc_op.wrap_synth(alg_loop, n);

                        %then put them back
                        %this is because the reduced-order-control
                        %implementation is done with a limited set of
                        %available channels
                        alg_psi.nw = alg_psi.nw + nwp_orig;
                        alg_psi.nz = alg_psi.nz + nzp_orig;
                        alg_psi.nwp = alg_psi.nwp - nwp_orig;
                        alg_psi.nzp = alg_psi.nzp - nzp_orig;
                    end
                    
                end
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
            %type of dynamical system (e.g. LTI, switched)
            tp = obj.type;
        end

        function ds = get_discount(obj)
            %should only certain modes be discounted
            ds = obj.discount;
        end
        

        %% getters
        
        %must define get_P, get_K

        function [Aa, B1, B2, C1, D11, D12, C2, D21, D22] = ss_zy_wu(obj, param)
            %get state space matrices at the current parameter values
            
            [Aa, B1, B2, C1, D11, D12, C2, D21, D22] = obj.P.ss_zy_wu();
        end

        function sys_alg = get_alg(obj, param)
            %close the loop of the algorithm
            %
            %Args: 
            %   param: structure of parameters
            %Return:
            %   sys_alg: the dynamical system interfacing the operators
            if nargin < 2
                param = [];
            end
            Pcurr = obj.get_P(param);

            if isnumeric(param) && isempty(param)
                param = struct('mode', 1);
            end
            
            if isempty(obj.K)
                sys_alg = Pcurr;
            else

                Kcurr = obj.get_K(param);
                sys_alg = lft(Pcurr, Kcurr);
            end
        end

        function op_out = get_op(obj, index)
            %get the operator at index
            %Args:
            %   index: the index

            op_out = obj.op{obj.bind(index)};
        end
        
        function pow = discount_schedule(obj, ordermax)
            %DISCOUNT_SCHEDULE exponential weights encountered when
            %applying the FIR filters
            %
            %Args:
            %   ordermax: maximum order of the IQCs
            %
            %Return:
            %   pow: Exponent sequence of discounts
            %
            % Example: 
            %   [0; 1 ; 2] -> rho.^[0; 1; 2] for uniform exponential stability
            
            %This becomes relevant when performing shuffled systems
            %(override on switched systems) 
            %
            %TODO: switched systems
            

            if obj.discount
                pow = -(0:ordermax)';
            end


        end


        %% for simulation

        function [y, u] = get_internal_signals(obj, param, x_all, w_all)
            %extract the internal signals from the interconnection (y, u) 
            %using the well-posedness expression
            %
            %Args: 
            %   param: structure of parameters
            %   x_all:      all states of network and controller
            %   w_all:      all inputs to the network (except u)            
            %Return:
            %   y:  input to controller/output of plant
            %   u:  output of controller/intput to plant
            

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


        function [Sbeta, Rbeta] = get_tracked_opt(obj, param)
            %GET_TRACKED_OPT get the tracked position of the optimal
            %solution
            %
            % :math:`\eta^*_{k+1} = S_\beta \eta^*, \beta^*_{k} = R_\beta
            % \eta_k`.
            %
            %
            % Args:
            %   param: structure of parameters
            %   
            %Returns:
            %   Sbeta: exosystem for optimal solution
            %   Rbeta: output of optimal solution
            
            if isempty(obj.tracking)
                Sbeta = 1;
                Rbeta = 1;
            else
                Sbeta = obj.tracking.Sbeta;
                Rbeta = obj.tracking.Rbeta;
            end

            c = obj.op{1}.c;
            Sbeta = kron(Sbeta, eye(c));
            Rbeta = kron(Rbeta, eye(c));
        end

        
        function mode_next = next_mode(obj, mode)
            %next mode in switching
            
            mode_next = 1;


        end

        function N = get_consensus_weighted(obj, op, bind)
            %GET_CONSENSUS_WEIGHTED create the consensus matrix
            %weight by the number of times the operator appears in bind
            %
            %Args:
            %   op: cell of operator
            %   bind: repeated patterns
            %Returns:
            %   N: consensus matrix
            if nargin < 2
                op = obj.op;
                bind = obj.bind;
            end

            N = get_consensus(obj, op, bind);

            [gc, gt] = groupcounts(reshape(bind, [], 1));

            weights = 1./(gc(bind));

            Nw = diag(weights) * N;


        end

        function N = get_consensus(obj, op, bind)
            %GET_CONSENSUS create the consensus matrix
            %for the regulation condition
            %
            %Args:
            %   op: cell of operator
            %   bind: repeated patterns
            %Returns:
            %   N: consensus matrix
            if nargin < 2
                op = obj.op;
                bind = obj.bind;
            end

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

            % c = op{1}.c;


            %index based on the bind 
            nbind = length(bind);
            Bind = full(sparse(1:nbind, bind, ~EQ(bind), nbind, nop));


            N = Bind * N0;
            % N = kron(N, eye(c));
        end    


        % function [obj, iwp, izp, supply] = add_ergodic_cert(obj, c)
        %     %ADD_ERGODIC_CERT certificate of ergodic convergence (function
        %     %value suboptimality). Used in conjunction with the op_sml.ERGODIC 
        %     % Section 4.1.2 (eq (32)) of https://arxiv.org/pdf/2302.06713
        %     %Args:
        %     %   c: kronecker lift of coordinate s            
        %     %Returns:
        %     %   iwp: new performance input indices
        %     %   izp: new performance output indices
        %     %   supply: IQZ
        % 
        %     [u, indbind] = unique(bind);
        %     nop = length(bind);
        %     nopu = length(obj.op);
        % 
        %     Nu = obj.get_consensus(obj.op, 1:nopu);
        % 
        %     ind_w = indbind;
        %     [obj.P, iwp] = obj.P.add_oracle_input(obj, indbind, []);
        % 
        % end

        function [iqc_curr, vars_curr,cons_curr] = create_iqc(obj, index, cons, order);
            %CREATE_IQC form the iqc for the current operator in the
            %system description
            %
            %Args:
            %   index (int): index of the operator
            %   cons: accumulated constraints
            %   order:  order of the operator: [causal order, noncausal order],
            %   or scalar for causal order
            %   reps: number of repeated evaluations (in bind)
            %Returns:
            %   iqc:    a valid iqc for the operator
            %   vars:   variables of the problem
            %   cons:   constraints in the problem (in terms of the
            %           variables directly)

            rep_curr = nnz(obj.bind == index);

            
            [iqc_curr, vars_curr,cons_curr] = obj.op{index}.create_iqc(cons, order{index}, rep_curr);
        end

        function sys_sim = export_sim(obj, op_sim)
            % export the system for use in simulation
            % with the operators (for iqcs) replaced by operators (in
            % op_sim)
            %
            % Args:
            %   op_sim: operators for simulation
            % Return:
            %   sys_sim: system for use in alg_sim
            sys_sim = obj;
            sys_sim.op = op_sim;

            %carry over the coordinate dimensions
            for i = 1:length(obj.op)
                sys_sim.op{i}.c = obj.op{i}.c;
            end


        end


    end

    methods (Abstract)
        get_K(obj, param) %get the controller at the current parameter values
        get_P(obj, param) %get the network at the current parameter values

        nxn(obj)    %number of network states
        nxi(obj)    %number of controller states


        build_plant() %form the plant to certify dissipation inequalities
    
    end
end

