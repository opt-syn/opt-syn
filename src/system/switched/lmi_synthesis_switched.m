classdef lmi_synthesis_switched < lmi_synthesis_interface
    %LMI_SYNTHESIS_SWITCHED synthesis LMIs for algorithmic interconnections
    %involving switched linear networks and controllers
    %
    %
    %examples include time-varying delays or coordinate updates
    %
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(mode(k))    Bw(mode(k))    Bwp(mode(k))   Bu(mode(k))  ][x(k)]   state transition
    % [z(k)  ] = [Cz(mode(k))   Dzw(mode(k))   Dzwp(mode(k))  Dzu(mode(k)) ][w(k)]   output to oracle
    % [zp(k) ] = [Czp(mode(k))  Dzpw(mode(k))  Dzpwp(mode(k)) Dzpu(mode(k))][wp(k)]  output to performance
    % [y(k) ]  = [Cy(mode(k))   Dyw(mode(k))   Dywp(mode(k))  Dyu(mode(k)) ][u(k)]   output to controller
    %   
    %
    %   Implemented
    %       stability
    %       e2e
    %       quad
    %       p2p
    %
    %   TODO:
    %       h2      
    %       e2p
    %       
    %    

    methods
        function obj = lmi_synthesis_switched(sys,config)
            %LMI_SYNTHESIS_SWITCHED undefined
            %   undefined
            obj@lmi_synthesis_interface(sys, config);
        end

        %% definition of variables and helpers

        function ns = Nss(obj)
            %NSS: Number of subsystems            
            ns = obj.sys.Nss;
        end

        function cm = common(obj)
            cm = obj.config.switched.common;
        end

        function [vars_diss, cons]= create_vars_storage(obj, cons, alg_psi, name)
            %create_vars_storage create variables for the dissipation
            %constraints. One for each subsystem
            %
            %
            %a cell of G(s) functions
            if nargin < 4
                name = [];
            end


            

            vars_diss = struct;


            
            % GX_cell = cell(obj.Nss, 1);
            % GY_cell = cell(obj.Nss, 1);
            % S_cell = cell(obj.Nss, 1);

            

            if obj.common
                %common storage function among all subsystems
                [GX, GY, cons] = obj.define_storage_G(cons, alg_psi{1}, '');
                vars_diss.JX = GX;
                vars_diss.JY = GY;
                
                n = ssize(GX, 1);

                JS = eye(n);
    
                vars_diss.JS = JS;

                for i = 1:obj.Nss
                    GX_cell{i} = GX;
                    GY_cell{i} = GY;
                    GS_cell{i} = eye(n);
                end

            else

                %common (slack) function
                [JX, JY, cons] = obj.define_storage_G(cons, alg_psi{1}, '_slack');
                vars_diss.JX = JX;
                vars_diss.JY = JY;

                n = ssize(JX, 1);
                vars_diss.JS = lmim(['JS_slack', name], n, n, 'full');                

                
                %per-mode storage functions, coupled by the mode-independent slack
                GX_cell = cell(obj.Nss, 1);
                GY_cell = cell(obj.Nss, 1);
                GS_cell = cell(obj.Nss, 1);
                              


                for i = 1:obj.Nss
                    [GX, GY, cons] = obj.define_storage_G(cons, alg_psi{i}, num2str(i));
                    n = ssize(GX, 1);
                    GX_cell{i} = GX;
                    GY_cell{i} = GY;
                    GS_cell{i} = lmim(['GS_slack', name], n, n, 'full');
                end

            end

            
            vars_diss.GX = GX_cell;
            vars_diss.GY = GY_cell;
            vars_diss.GS  = GS_cell;

        end


        function [vars_K, cons] = create_vars_controller(obj, cons, alg_psi, name)
            %CREATE_VARS_CONTROLLER create the nonlinearly-transformed
            %controller matrices

            %get the dimensions

            vars_K = cell(obj.Nss, 1);

            if nargin < 4
                name = [];
            end

            D_mask = obj.get_D_mask();            
            for i = 1:obj.Nss
                name_curr = [name, '_', num2str(i)];
                alg_curr = alg_psi{i};
                D_mask_curr = D_mask{i};
                [vars_K{i}, cons] = create_vars_controller@lmi_synthesis_interface(obj, cons, alg_curr, name_curr, D_mask_curr);
            end
            
        end

        function D_mask = get_D_mask(obj)
            %GET_D_MASK get the direct feedthrough terms

            %the sparsity-constrained term for internal model control            
            D_mask_0 = obj.config.syn.D_mask;

            

            D_mask_default  = tril(ones(length(obj.sys.bind)));
            % Handle the case when D_mask_0 is empty


            if ~iscell(D_mask_0)
                D_mask_0_orig = D_mask_0;
                D_mask_0 = cell(obj.Nss, 1);
                for i = 1:obj.Nss
                    if isempty(D_mask_0_orig)
                        D_mask_0{i} = D_mask_default;
                    else
                        D_mask_0{i} = D_mask_0_orig;
                    end
                end
            end
            

            %WARNING: do a better conversion on the coordinate lifts
            c = obj.sys.op{1}.c;
            D_mask = cell(obj.Nss, 1);
            for i = 1:obj.Nss
                D_mask{i} = kron(D_mask_0{i}, ones(c));
            end

        end

        function vars_inv= get_vars_involved(obj, vars, ind)
            %GET_VARS_INVOLVED get variables involved in the current mode

            vars_inv= struct;
            if ind==0
                vars_inv.diss.GX = vars.diss.JX;
                vars_inv.diss.GY = vars.diss.JY;
                vars_inv.diss.GS  = vars.diss.JS;
            else
                vars_inv.diss.GX = vars.diss.GX{ind};
                vars_inv.diss.GY = vars.diss.GY{ind};
                vars_inv.diss.GS  = vars.diss.GS{ind};
            end



            if ind 
                vars_inv.reg.Pi = vars.reg.Pi{ind};
                vars_inv.reg.Gam = vars.reg.Gam{ind};
                vars_inv.reg.Phi = vars.reg.Phi{ind};
            else
                vars_inv.reg =[];
            end

        end

        %% main call
        function [vars, cons, objective, con_M] = cons_dynamic(obj, vars, cons, diss)
            %CONS form the dissipation and sign constraints
            %
            %Input:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss:   structure describing the problem
            %       plant:  system to control
            %       spec:   performance specification           
            %       target: whether the performance measure should be optimized
            %               true:  soft constraint (e.g. Schur complement
            %                                       formulation)
            %               false: hard constraint            
            %       ind_curr:  the index of the current subsystem
            %       ind_next:  the index of the next subsystem
            %
            %Output:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            


            %need to look up the right constraint            

            %Upper-levels: iterate over the systems
            objective = 0;

            ds = obj.sys.get_discount();
            if obj.common

                for i = 1:obj.Nss
                    %extract the information of subsystem i
                    diss_curr = diss;
                    diss_curr.plant = diss.plant{i};
                    diss_curr.ind_curr = (i);
                    diss_curr.ind_next = (i);

                    if ds(i)
                        diss_curr.rho = diss.rho;
                    else
                        diss_curr.rho = 1;
                    end


                    [cons, objective_curr, con_M] = obj.con_dynamic_single(vars, cons, diss_curr);


                    %TODO: take the max over the different subsystems
                    %but the same objective is sent to each subsystem, so it's
                    %all the same? Check this
                    if i==1
                        objective = objective + objective_curr;
                    end

                end 

            else
                [src, dst] = obj.sys.get_arcs();
                Narcs = length(src);

                for i = 1:Narcs
                    %extract the information of subsystem i
                    diss_curr = diss;
                    diss_curr.plant = diss.plant{src(i)};
                    diss_curr.ind_curr = src(i);
                    diss_curr.ind_next = dst(i);

                    
    
                    if ds(src(i))
                        diss_curr.rho = diss.rho;
                    else
                        diss_curr.rho = 1;
                    end
    
                    [cons, objective_curr, con_M] = obj.con_dynamic_single(vars, cons, diss_curr);
    
    
                    %TODO: take the max over the different subsystems
                    %but the same objective is sent to each subsystem, so it's
                    %all the same? Check this
                    if i==1
                        objective = objective + objective_curr;
                    end
    
                end 

                %impose sign constraint    
                for i = 1:obj.Nss
                    vcurr = obj.get_vars_involved(vars, i);
                    Gcurr = obj.get_storage(vcurr.diss, vcurr.reg);
                    cons = obj.con_terminal(Gcurr, cons, [], diss.iqc_rob);
                end
            end

        end

        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance



            %get the variables of the problem
            vslack = obj.get_vars_involved(vars, 0);            
            vcurr = obj.get_vars_involved(vars, diss.ind_curr);
            vnext = obj.get_vars_involved(vars, diss.ind_next);
            

            Gslack = obj.get_storage(vslack.diss, vslack.reg);
            Gcurr = obj.get_storage(vcurr.diss, vcurr.reg);
            Gnext = obj.get_storage(vnext.diss, vnext.reg);



            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            P = obj.reg.connect_model(diss.plant, diss.ind_curr, diss.rho);            

            sys_cl = obj.system_closed_loop(P, vslack.diss, vslack.reg, vars.K{diss.ind_curr});
            
            %index the quadratic specification
            np = diss.iqc_rob.np;
            nq = diss.iqc_rob.nq;

            M_quad_rob = quad_objective_decomp(diss.iqc_rob.M, 1:np, np + (1:nq));
            [M_quad_spec, objective] = diss.spec.supply_quad(vars.spec);

            quad = obj.merge_quad(M_quad_rob, M_quad_spec);
           
            %formulation from ParDynSyn notes (parametric dynamic
            %synthesis)


            %the quadratic objective
            supp_b = obj.supply_block(sys_cl, quad);

            %the storage
            Gcurr_slack = Gslack + Gslack' - Gcurr;
            stor_b = obj.storage_block(sys_cl, quad, Gcurr_slack, Gnext);

            %the dynamics
            dyn_b = obj.dynamics_block(sys_cl, quad);
            
            %wrap it all together
%             objective = 0;

            con_M = -(stor_b + supp_b + dyn_b);


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 
            
        end

        function cons = con_spread(obj, cons, vars)
            %CON_SPREAD increase numerical conditioning by separating the 
            %primal and dual blocks
            %invoke this over multiple subsystems
            if ~obj.config.syn.reduced_order
                for i = 1:obj.Nss
                    cons = obj.con_spread_single(cons, vars.diss.GX{i}, vars.diss.GY{i});
                end
            end
        end

        %TODO: e2e_target

        %% Recovery
        % Call the recovery method to finalize the synthesis process
        % [cons, objective, con_M] = obj.recovery(vars, cons, diss);

        function [sol] = recover_subcontroller(obj, alg_psi, P_trans, sol)
            %RECOVER_SUBCONTROLLER recover the subcontroller of the current
            %mode/control
            %
            %
            %Input:
            %
            %Output:
            %   K_feed: the subcontroller with direct feedthrough, before
            %           exponential discounting    
            %(not yet exponentially undiscounted, this happens later)


            vars_rec = sol.vars;
            rho = sol.rho;
                    
            %recover the controller
            [K_nofeed] = recover_subcontroller_warp(obj, P_trans, vars_rec);

            %package it up
            K_report = cell(obj.Nss, 1);

            ds = obj.sys.get_discount;

            for i = 1:obj.Nss

                if ds(i)
                    curr_rho = rho;
                else
                    curr_rho = 1;
                end
                
                model = obj.reg.get_model(i, vars_rec.reg);
    
                K_report = obj.K_alg_report(P_trans{i}, K_nofeed{i}, model, curr_rho);

                %form the algorithm
                sol.alg_trans{i} = K_report.alg_trans;
                sol.alg{i} = lft(obj.sys.P{i}, K_report.K);
                sol.model{i} = K_report.model;           
                sol.K{i}= K_report.K;
                sol.K_sub{i} = K_report.K_sub;                

            end


            sol.gain = obj.validate_recovery_gain(sol.alg_trans, sol.iqc_op_all);
        end

        function [K_nofeed] = recover_subcontroller_warp(obj, P_trans, vars_rec)

            %RECOVER_SUBCONTROLLER_WARP recover the nonlinearly warped
            %controller 
            %dynamics and indexers


            %for debugging
            % G = obj.get_storage(sol.vars.diss, sol.vars.reg);

            %this is the (nonlinearly-warped) system that is certified as
            %possessing the desired performance and robustness
            %specifications
            % sys_cl = obj.system_closed_loop(P_trans, sol.vars.diss, sol.vars.reg, sol.vars.K);

            % sys_cal = ss(G \ Acl, G \ Bcl, Ccl, Dcl, 1);


            %extract the storage variables and the factorizations       

            K_nofeed = cell(obj.Nss, 1);        



            JS = vars_rec.diss.JS;
            JX = vars_rec.diss.JX;
            JY = vars_rec.diss.JY;

            J = JS - JX * JY;
            [Up, Sig, Vp] = svd(J);

            % U = Up*Sig;
            ssig = sqrt(Sig);
            srsig = diag(1./(diag(ssig)));


            V = Vp*ssig;
            U = Up*ssig;

            Uinv = srsig*Up';
            Vinv = srsig*Vp';
            
            
            %get the indexers
            Pt = P_trans{1};            

            iz = [Pt.index_z(), Pt.index_zp()];
            iw = [Pt.index_w(), Pt.index_wp()];
            iu = Pt.index_u();
            iy = Pt.index_y();           

            nz = length(iz);
            nw = length(iw);
            nu = length(iu);
            ny = length(iy);

            for i = 1:obj.Nss
                [A, B, C, D] = ssdata(P_trans{i});
                Ak = vars_rec.K{i}.A;
                Bk = vars_rec.K{i}.B;
                Ck = vars_rec.K{i}.C;
                Dk = vars_rec.K{i}.D;
    
                
                n = ssize(Ak, 1);                
    
    
                %controller recovery
    
                Lblock = [Uinv, -Uinv*JX*B(:, iu);
                    zeros(nu, size(V, 2)), eye(nu)];
    
                
                Cblock = [Ak - JX*A*JY, Bk;
                    Ck, Dk];
  
    
                Rblock = [Vinv', zeros(size(V, 2), ny);
                    -C(iy, :)*JY*Vinv', eye(ny)];
    
    
                % Kblock0 = inv(LblockI)* (Cblock) * inv(RblockI);
                % Kblock1 = LblockI) \ Cblock) * inv(RblockI);
                % Kblockinv = RblockI') \ LblockI) \ Cblock)')';
    
                Kblock = Lblock * Cblock * Rblock;
    
                %extraction and exponential weighting
                Ac = Kblock(1:n, 1:n);
                Bc = Kblock(1:n, n+1:end);
                Cc = Kblock(n+1:end, 1:n);
                Dc = Kblock(n+1:end, n+1:end);
    
                K_nofeed{i}= ss(Ac, Bc, Cc, Dc, 1);                
            end
        end

        function gain = validate_recovery_gain(obj, alg_trans, iqc_op_all)
            %VALIDATE_RECOVERY validate that the system obeys the stability
            %constraint (TODO: performance specs)

            %not yet supported
            gain = 0;
        end
    end
end