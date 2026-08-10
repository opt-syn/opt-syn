classdef lmi_synthesis_periodic < lmi_synthesis_switched
    %LMI_SYNTHESIS_PERIODIC synthesis LMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
    
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(k)    Bw(k)    Bwp(k)   Bu(k)  ][x(k)]   state transition
    % [z(k)  ] = [Cz(k)   Dzw(k)   Dzwp(k)  Dzu(k) ][w(k)]   output to oracle
    % [zp(k) ] = [Czp(k)  Dzpw(k)  Dzpwp(k) Dzpu(k)][wp(k)]  output to performance
    % [y(k)  ] = [Cy(k)   Dyw(k)   Dywp(k)  Dyu(k) ][u(k)]   output to controller
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

    methods
        function obj = lmi_synthesis_periodic(sys,config)
            %LMI_SYNTHESIS_PERIODIC constructor
            obj@lmi_synthesis_switched(sys, config);
        end

        %% definition of variables and helpers


        function [vars_diss, cons]= create_vars_storage(obj, cons, alg_psi, name)
            %create_vars_storage create variables for the dissipation
            %constraints. One for each subsystem
            %Args:                   
            %   cons:   accumulated constraints
            %   specs: performance specifications
            %
            %Returns:            
            %  vars_spec:   variables for performance specification
            %   cons:   accumulated constraints
            

            if nargin < 4
                name = [];
            end           


            GX_cell = cell(obj.Nss, 1);
            GY_cell = cell(obj.Nss, 1);
            GS_cell = cell(obj.Nss, 1);

            if obj.common
                %common storage function among all subsystems
                
                [GX, GY, cons] = obj.define_storage_G(cons, alg_psi{1}, '');
                n = ssize(GX, 1);
                

                G = vars_diss.G;
                GX_cell = cell(obj.Nss, 1);
                for i = 1:obj.Nss
                    GX_cell{i} = GX;
                    GY_cell{i} = GY;
                    GS_cell{i} = eye(n);
                end

            else
                %define a storage function for each subsystem

                for i = 1:obj.Nss
                    [GX, GY, cons] = obj.define_storage_G(cons, alg_psi{i}, num2str(i));
                    n = ssize(GX, 1);
                    GX_cell{i} = GX;
                    GY_cell{i} = GY;
                    GS_cell{i} = eye(n);
                end

            end

            vars_diss = struct;
            vars_diss.GX = GX_cell;
            vars_diss.GY = GY_cell;
            vars_diss.GS  = GS_cell;

        end

        function vars_inv= get_vars_involved(obj, vars, ind)
            %GET_VARS_INVOLVED get variables involved in the current mode
            %Args:
            %   vars:   variables of the problem        
            %   ind:    index of subsystem/mode
            %Returns:
            %   vars_inv:     variables (diss, reg) at subsystem ind


            vars_inv= struct;
            vars_inv.diss.GX = vars.diss.GX{ind};
            vars_inv.diss.GY = vars.diss.GY{ind};
            vars_inv.diss.GS  = vars.diss.GS{ind};

            vars_inv.reg.Pi = vars.reg.Pi{ind};
            vars_inv.reg.Gam = vars.reg.Gam{ind};
            vars_inv.reg.Phi = vars.reg.Phi{ind};

        end

        %% main call
        function [vars, cons, objective] = cons_dynamic(obj, vars, cons, diss)
            %CONS form the dissipation and sign constraints
            %
            %Args:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss (diss_data):   structure describing the dissipation
            %       constraint
            %
            %Returns:
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

        %% performance specifications

        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance
            %
            %Args:                   
            %   cons:   accumulated constraints
            %   specs: performance specifications
            %
            %Returns:            
            %  vars_spec:   variables for performance specification
            %   cons:   accumulated constraints
                


            %get the variables of the problem
            vcurr = obj.get_vars_involved(vars, diss.ind_curr);
            vnext = obj.get_vars_involved(vars, diss.ind_next);
            
            Gcurr = obj.get_storage(vcurr.diss, vcurr.reg);
            Gnext = obj.get_storage(vnext.diss, vnext.reg);



            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            rhou = obj.used_rho(diss);
            P = obj.connect_model(diss, rhou);            

            

            %dynamics constraint
            vars_diss = vcurr.diss;
            vars_diss.GX = vnext.diss.GX;
            
            sys_cl = obj.system_closed_loop(P, vars_diss, vars.reg, vars.K{diss.ind_curr}, diss.rho);
            
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
            stor_b = obj.storage_block(sys_cl, quad, Gcurr, Gnext);

            %the dynamics
            dyn_b = obj.dynamics_block(sys_cl, quad);
           

            con_M = -(stor_b + supp_b + dyn_b);


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint            
            cons = obj.con_terminal(Gcurr, cons, [], diss.iqc_rob);
        end

        %% helper functions       
        

 

        %% Recovery
        % Call the recovery method to finalize the synthesis process

        function [sol] = recover_subcontroller(obj, alg_psi, P_trans, sol)
            %RECOVER_SUBCONTROLLER recover the subcontroller of the current
            %mode/control
            %
            %
            %Args:
            %   alg_psi:   the filtered algorithmic interconnection
            %   P_trans:    the transformed generalized plant before IQC
            %   sol: solution structure
            %
            %Returns:
            %   sol: solution structure
            

            vars_rec = sol.vars;

        
            %recover the controller
            if obj.config.gen.same_rho
                rho_common = sol.rho;
            else
                rho_common = 1;
            end
            [K_nofeed, Gcl, Ycl] = recover_subcontroller_warp(obj, P_trans, vars_rec);


            %package it up
            K_report = cell(obj.Nss, 1);
            for i = 1:obj.Nss

                model = obj.reg.get_model(i, vars_rec.reg);
    
                K_report = obj.K_alg_report(P_trans{i}, K_nofeed{i}, model, rho_common);

                %form the algorithm
                sol.cert.alg_trans{i} = K_report.alg_trans;
                sol.cert.alg{i} = lft(obj.sys.P{i}, K_report.K);
                sol.cert.model{i} = K_report.model;           
                sol.cert.K{i}= K_report.K;
                sol.cert.K_sub{i} = K_report.K_sub;                

            end
            sol.cert.Gcl = Gcl;
            sol.cert.Ycl = Ycl;

            n = sol.cert.alg_trans{1}.dump_dim();
            sol.cert.alg_trans = genplant_poly(sol.cert.alg_trans, n);

            

            sol.gain = obj.validate_recovery_gain(sol.cert.alg_trans, sol.cert.iqc_op_all);
        end

        function [K_nofeed, Gcl, Ycl] = recover_subcontroller_warp(obj, P_trans, vars_rec)
            %RECOVER_SUBCONTROLLER_WARP recover the nonlinearly warped
            %controller 
            %dynamics and indexers
            %Args:
            %   alg_psi:   the filtered algorithmic interconnection
            %   P_trans:    the transformed generalized plant before IQC
            %   sol: solution structure
            %
            %Output:
            %   K_nofeed: subcontroller without direct feedthrough
            %   Gcl:    closed-loop storage matrix (original)
            %   Ycl:    similarity transformation/nonlinear warping


            

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

            Gcl = cell(obj.Nss, 1);
            Ycl = cell(obj.Nss, 1);          

            K_nofeed = cell(obj.Nss, 1);        

            % G = obj.get_storage(sol.vars.diss, sol.vars.reg);

            for i = 1:obj.Nss
                Y{i} = vars_rec.diss.GY{i};
                X{i} = vars_rec.diss.GX{i};

                S = vars_rec.diss.GS{i};

                J = S - X{i} * Y{i};
                [Up, Sig, Vp] = svd(J);

                n = size(Y{i}, 1);
                
                Ycl{i} = [Y{i}, eye(n); Vp', zeros(n)];

                iYcl = inv(Ycl{i}); 
                Gcl{i} = iYcl' * [Y{i}, eye(n); eye(n), X{i}] * iYcl;
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

        function gain = validate_recovery_gain(obj, alg_trans, iqc_op_all)
            %VALIDATE_RECOVERY validate that the system obeys the stability
            %Args:
            %   alg_trans: the plant with confirmed performance by LMIs
            %   iqc_op_all: all IQCs
            %Return:
            %   gain:   [Passivity index, H-infinity index].


            %use the monodromy system to get specs
            % n = alg_trans{1}.dump_dim();
            alg_trans_lti = periodic_lift(alg_trans.P);
            c = obj.sys.op{1}.c;


            P = alg_trans_lti;


            M = kron(iqc_op_all.iqc.M, eye(c));
            M = (M + M')/2;
            nw = floor(size(M, 1)/2);

            M11 = M(1:nw, 1:nw);
            M12 = M(nw + (1:nw), 1:nw);
            M22 = M(nw + (1:nw), nw + (1:nw));
            %is the constraint passive?
            is_passive = (norm(M11) + norm(M22) + norm(M12 - eye(nw)))==0;
            is_hinf = (norm(M11-eye(nw)) + norm(M22+eye(nw)) + norm(M12))==0;

            nw_lift = size(P.D, 1);       
            E=eye(nw_lift);

            if is_passive
                gain_passive = -getPassiveIndex(-P, 'input');

                
                
                Tinf=[E sqrt(2)*E;sqrt(2)*E E];
                P_inf = lft(Tinf,P,nw_lift,nw_lift);

                gain_inf = norm(P_inf, 'inf');
            elseif is_hinf
                gain_inf = norm(P, 'inf');

                
                Tpass = [-E sqrt(2)*E;sqrt(2)*E -E];
                Ppass = lft(Tpass,P,nw_lift,nw_lift);

                gain_passive = -getPassiveIndex(-Ppass, 'input');
            else
                %TODO: advanced validation
                warning('Customized validation is not yet implemented')
                gain_inf = 0;
                gain_passive = 0;
            end

            gain = [gain_passive, gain_inf];            

        end

    end
end