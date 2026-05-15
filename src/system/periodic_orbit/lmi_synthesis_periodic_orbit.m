classdef lmi_synthesis_periodic_orbit < lmi_synthesis_periodic
    %LMI_SYNTHESIS_PERIODIC synthesisLMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(k)    Bw(k)    Bwp(k)   Bu(k)  ][x(k)]   state transition
    % [z(k)  ] = [Cz(k)   Dzw(k)   Dzwp(k)  Dzu(k) ][w(k)]   output to oracle
    % [zp(k) ] = [Czp(k)  Dzpw(k)  Dzpwp(k) Dzpu(k)][wp(k)]  output to performance
    % [zp(k) ] = [Cy(k)   Dyw(k)   Dywp(k)  Dyu(k) ][u(k)]   output to controller
    %
    %A(k) = A(k+T) for some known time T
    %
    %instances of these algorithms include cyclic coordinate descent
    %methods. Periodic systems can also be unrolled into an LTI system
    %(monodromy methods): a single large LMI system rather than multiple 
    % coupled smaller LMI systems
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

    properties
        opts = struct("COMMON", false);
    end

    methods
        function obj = lmi_synthesis_periodic(sys,config)
            %LMI_SYNTHESIS_PERIODIC undefined
            %   undefined
            obj@lmi_synthesis_interface(sys, config);
        end

        %% definition of variables and helpers

        function ns = Nss(obj)
            %NSS: Number of subsystems            
            ns = obj.sys.Nss;
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


            


            GX_cell = cell(obj.Nss, 1);
            GY_cell = cell(obj.Nss, 1);
            S_cell = cell(obj.Nss, 1);

            if obj.opts.COMMON
                %common storage function among all subsystems
                
                [GX, GY, cons] = obj.define_storage_G(cons, alg_psi{1}, '');
                n = ssize(GX, 1);
                % vars_diss= struct('GX', GX, 'GY', GY, 'S', eye(n));

                G = vars_diss.G;
                GX_cell = cell(obj.Nss, 1);
                for i = 1:obj.Nss
                    GX_cell{i} = GX;
                    GY_cell{i} = GY;
                    S_cell{i} = eye(n);
                end

            else
                %define a storage function for each subsystem

                for i = 1:obj.Nss
                    [GX, GY, cons] = obj.define_storage_G(cons, alg_psi{i}, num2str(i));
                    n = ssize(GX, 1);
                    GX_cell{i} = GX;
                    GY_cell{i} = GY;
                    S_cell{i} = eye(n);
                end

            end

            vars_diss = struct;
            vars_diss.GX = GX_cell;
            vars_diss.GY = GY_cell;
            vars_diss.S  = S_cell;

        end


        function [vars_K, cons] = create_vars_controller(obj, cons, alg_psi, name)
            %CREATE_VARS_CONTROLLER create the nonlinearly-transformed
            %controller matrices

            %get the dimensions

            vars_K = cell(obj.Nss, 1);

            if nargin < 4
                name = [];
            end

            for i = 1:obj.Nss
                name_curr = [name, '_', num2str(i)];
                alg_curr = alg_psi{i};
                [vars_K{i}, cons] = create_vars_controller@lmi_synthesis_interface(obj, cons, alg_curr, name_curr);
            end
            
        end

        function vars_inv= get_vars_involved(obj, vars, ind)
            %GET_VARS_INVOLVED get variables involved in the current mode

            vars_inv= struct;
            vars_inv.diss.GX = vars.diss.GX{ind};
            vars_inv.diss.GY = vars.diss.GY{ind};
            vars_inv.diss.S  = vars.diss.S{ind};

            vars_inv.reg.Pi = vars.reg.Pi{ind};
            vars_inv.reg.Gam = vars.reg.Gam{ind};
            vars_inv.reg.Phi = vars.reg.Phi{ind};

        end

        %% main call
        function [cons, objective, con_M] = cons_dynamic(obj, vars, cons, diss)
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
            for i = 1:obj.Nss
                %extract the information of subsystem i
                diss_curr = diss;
                diss_curr.plant = diss.plant{i};
                diss_curr.ind_curr = i;
                diss_curr.ind_next = 1+mod(i, obj.Nss);


                [cons, objective_curr, con_M] = obj.con_dynamic_single(vars, cons, diss_curr);


                %TODO: take the max over the different subsystems
                %but the same objective is sent to each subsystem, so it's
                %all the same? Check this
                if i==1
                    objective = objective + objective_curr;
                end
            end         

        end

        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance



            %get the variables of the problem
            vcurr = obj.get_vars_involved(vars, diss.ind_curr);
            vnext = obj.get_vars_involved(vars, diss.ind_next);
            
            Gcurr = obj.get_storage(vcurr.diss, vcurr.reg);
            Gnext = obj.get_storage(vnext.diss, vnext.reg);



            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            P = obj.reg.connect_model(diss.plant, diss.ind_curr, diss.rho);

            vars_diss = vcurr.diss;
            vars_diss.GX = vnext.diss.GX;

            sys_cl = obj.system_closed_loop(P, vars_diss, vars.reg, vars.K{diss.ind_curr});
            
            %index the quadratic specification
            vars_spec = vars.spec{diss.spec.id};
            M_quad = -obj.merge_spec_M(diss.iqc_rob, diss.spec, vars_spec);


            if isempty(diss.spec.izp)
                ind_p = 1:(diss.iqc_rob.nz);
                ind_q = diss.iqc_rob.nz + (1:(diss.iqc_rob.nw));
            else
                ind_p = 1:(diss.iqc_rob.nz + diss.spec.izp);
                ind_q = (diss.iqc_rob.nz + diss.spec.izp) + (1:(diss.iqc_rob.nw + diss.spec.iwp));
            end
            
            quad = obj.quad_objective(M_quad, ind_p, ind_q);
                       
            %formulation from ParDynSyn notes (parametric dynamic
            %synthesis)


            %the quadratic objective
            supp_b = obj.supply_block(sys_cl, quad);

            %the storage
            stor_b = obj.storage_block(sys_cl, quad, Gcurr, Gnext);

            %the dynamics
            dyn_b = obj.dynamics_block(sys_cl, quad);
            
            %wrap it all together
            objective = 0;

            con_M = stor_b + supp_b + dyn_b;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint            
            cons = obj.con_terminal(Gcurr, cons, [], diss.iqc_rob);
        end

        %TODO: e2e_target

        %% Recovery
        % Call the recovery method to finalize the synthesis process
        % [cons, objective, con_M] = obj.recovery(vars, cons, diss);

        function gain = validate_recovery_gain(obj, alg_trans, iqc_op_all);

            %Validation of the LMI
            %not performed at the moment.
            gain = 0;

        end

        function [sol] = recover_subcontroller(obj, P_trans, sol)
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
            for i = 1:obj.Nss

                model = obj.reg.get_model(i, vars_rec.reg);
    
                K_report = obj.K_alg_report(P_trans{i}, K_nofeed{i}, model, rho);

                %form the algorithm
                sol.alg_trans{i} = K_report.alg_trans;
                sol.alg{i} = lft(obj.sys.P{i}, K_report.K);
                sol.model{i} = K_report.model;           
                sol.K{i}= K_report.K;
                sol.K_sub{i} = K_report.K_sub;                

            end


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
            Y = cell(obj.Nss, 1);
            X = cell(obj.Nss, 1);
            U = cell(obj.Nss, 1);
            V = cell(obj.Nss, 1);
            Uinv = cell(obj.Nss, 1);
            Vinv = cell(obj.Nss, 1);            

            K_nofeed = cell(obj.Nss, 1);        

            for i = 1:obj.Nss
                Y{i} = vars_rec.diss.GY{i};
                X{i} = vars_rec.diss.GX{i};

                S = vars_rec.diss.S{i};

                J = S - X{i} * Y{i};
                [Up, Sig, Vp] = svd(J);

                % U = Up*Sig;
                ssig = sqrt(Sig);
                srsig = diag(1./(diag(ssig)));


                V{i} = Vp*ssig;
                U{i} = Up*ssig;

                Uinv{i} = srsig*Up';
                Vinv{i} = srsig*Vp';
            end
            
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
                inext = mod(i, obj.Nss) + 1;
    
    
                %controller recovery
    
                Lblock = [Uinv{inext}, -Uinv{inext}*X{inext}*B(:, iu);
                    zeros(nu, size(V{inext}, 2)), eye(nu)];
    
                
                Cblock = [Ak - X{inext}*A*Y{i}, Bk;
                    Ck, Dk];
  
    
                Rblock = [Vinv{i}', zeros(size(V{i}, 2), ny);
                    -C(iy, :)*Y{i}*Vinv{i}', eye(ny)];
    
    
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
    end
end