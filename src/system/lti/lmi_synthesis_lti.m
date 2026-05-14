classdef lmi_synthesis_lti < lmi_synthesis_interface
    %LMI_SYNTHESIS_LTI synthesis LMIs for algorithmic interconnections
    %involving linear-time-invariant (LTI) networks and controllers
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A    B     Bp     Bu  ][x(k)]   state transition
    % [z(k)  ] = [Cz   Dzw   Dzwp   Dzu ][w(k)]   input to oracle
    % [zp(k) ] = [Czp  Dzpw  Dzpwp  Dzpu][wp(k)]  performance specification
    % [y(k)  ] = [Cy   Dyw   Dywp   Dyu ][u(k)]   input to controller
    %
    % performance specification: wp -> zp from (spec)
    %
    %   Implemented:
    %
    %   TODO:
    %       stability
    %       e2e
    %       quad
    %       p2p
    %       h2      
    %       e2p
    %       
    %
    
    methods
        function obj = lmi_synthesis_lti(sys, config)
            %LMI_SYNTHESIS_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_synthesis_interface(sys, config);
        end       
        
        
        %% definition of variables and helpers
        function [vars_diss, cons]= create_vars_storage(obj, cons, alg_psi, name)
            %create_vars_storage create variables for the dissipation
            %constraints
            %
            %Input:
            %   cons:       accumulated constraints
            %   alg_psi:    the filtered algorithmic interconnection
            %   name:       a name for the variable

            if nargin < 4
                name = [];
            end

            [GX, GY, cons] = obj.define_storage_G(cons, alg_psi, name);
            vars_diss= struct('GX', GX, 'GY', GY);

        end

        %% stability (for testing)
        function [cons, objective, con_M] = stability_passive(obj, vars, cons, diss)
        % function [cons, objective, con_M] = stability(obj, vars, cons, diss)
            %certification of exponential stability
            
            G = obj.get_storage(vars.diss, vars.reg);

            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            P = obj.reg.connect_model(diss.plant, diss.rho);

            sys_cl = obj.system_closed_loop(P, vars.diss, vars.reg, vars.K);

            %only do this if the system is passive

            n = ssize(sys_cl.A, 1);
            nw = ssize(sys_cl.B, 2);
            nz = ssize(sys_cl.C, 1);
            
            % nt = length(ind_pos);

            dyn_block =  [G,  sys_cl.A, sys_cl.B;
            sys_cl.A', G, zeros(n, nz);
            sys_cl.B', zeros(nz, n), zeros(nz)];

            %supply block
            %       sp = [zeros(n), zeros(n), zeros(n, nz);
            % zeros(n), zeros(n), zeros(n, nz); %vars.ga * eye(n)
            % zeros(nz, n), Ccl, Dcl + eye(nz)*p.opts.pass_tol] * (-0.5);

            dissI = eye(nz)*obj.config.tol.input_diss;
            % dissI = zeros(nz);
              sp = [zeros(n), zeros(n), zeros(n, nz);
            zeros(n), zeros(n), zeros(n, nz); %vars.ga * eye(n)
            zeros(nz, n), sys_cl.C, sys_cl.D + dissI] * (-0.5);
            
            cost_block = sp + sp';
    
            con_M = dyn_block + cost_block;
    
            objective = 0;            


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint
            %change this up
            % cons = obj.con_terminal(G, cons, [], diss.iqc_rob);
        end


        %% Quadratic performance (infinite horizon)        
        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance
            
            %get the variables of the problem
            G = obj.get_storage(vars.diss, vars.reg);
            
            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            P = obj.reg.connect_model(diss.plant, diss.rho);

            sys_cl = obj.system_closed_loop(P, vars.diss, vars.reg, vars.K);
            
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
            stor_b = obj.storage_block(sys_cl, quad, G, G);

            %the dynamics
            dyn_b = obj.dynamics_block(sys_cl, quad);
            
            %wrap it all together
            objective = 0;

            con_M = stor_b + supp_b + dyn_b;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint            
            cons = obj.con_terminal(G, cons, [], diss.iqc_rob);
        end        



        function [cons, objective, con_M] = e2e_target(obj, vars, cons, diss)
            %E2E_TARGET: use a Schur complement to minimize the energy to
            %energy gain of the transfer function


            
                       %get the variables of the problem
            G = obj.get_storage(vars.diss, vars.reg);
            
            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            P = obj.reg.connect_model(diss.plant, diss.rho);

            sys_cl = obj.system_closed_loop(P, vars.diss, vars.reg, vars.K);
            
            %index the quadratic specification
            vars_spec = vars.spec{diss.spec.id};
            M_quad = -diss.iqc_rob.M;            


            ind_p = 1:(diss.iqc_rob.nz);
            ind_q = diss.iqc_rob.nz + (1:(diss.iqc_rob.nw));
            
            
            mu_l2 = vars.spec{diss.spec.id}.mu_l2;

            quad_rob = obj.quad_objective(M_quad, ind_p, ind_q);

            %adapt the quadratic objecitive for the e2e target
            nwp = length(diss.spec.iwp);
            nzp = length(diss.spec.izp);
            
            Q_e2e = -eye(nzp) * mu_l2;
            T_e2e = eye(nwp);
            S_e2e = zeros(nzp, nwp);
            U_e2e = -eye(nzp) * mu_l2;

            Q_new = blkdiag(quad_rob.Q, Q_e2e);
            T_new = blkdiag(quad_rob.T, T_e2e);
            S_new = blkdiag(quad_rob.S, S_e2e);
            U_new = blkdiag(quad_rob.U, U_e2e);

            quad = struct('Q', Q_new, 'T', T_new, 'S', S_new, 'U', U_new);

            
            %formulation from ParDynSyn notes (parametric dynamic
            %synthesis)


            %the quadratic objective
            supp_b = obj.supply_block(sys_cl, quad);

            %the storage
            stor_b = obj.storage_block(sys_cl, quad, G, G);

            %the dynamics
            dyn_b = obj.dynamics_block(sys_cl, quad);
            
            %wrap it all together
            objective = mu_l2;
            con_M = stor_b + supp_b + dyn_b;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint            
            cons = obj.con_terminal(G, cons, [], diss.iqc_rob);
            
            
           
        end

        %% Peak-to-Peak norm (at each finite horizon)

        function [cons, objective, con_M] = p2p(obj, vars, cons, diss)
            %p2p: certificate of peak to peak induced norm
            %
            % sup norm(zp, 2) / norm(wp, 2) <= objective

            % verification by Theorem 4 of https://www.sciencedirect.com/science/article/pii/S2405896323008194

            %TODO: fix exponential rate here
            %storage matrix

            error('LTI synthesis: p2p target not yet supported')
            G = vars.diss.G;
                      
            %terminal constraint
            X = diss.iqc_rob.X;
            nf = ssize(X);
            n = ssize(G, 1);
            Ef = [eye(nf); zeros(n-nf, nf)];

            X_f = Ef * X * Ef';

            

            %variables in the problem
            vars_spec = vars.spec{diss.spec.id};
            mu = vars_spec.mu_p2p;            

            %form the plant
            [plant_no_p, CDp] = obj.separate_performance_output(diss);


            %Block 1: without performance
            
            nwp = length(diss.spec.iwp);
            M_base = blkdiag(diss.iqc_rob.M, -mu * eye(nwp));
            
            sysb_1 = obj.sys_block(diss.plant, G, G);
            suppb_1 = obj.supply_block(plant_no_p, M_base);


            con_M_1 = sysb_1 + suppb_1;

            %Block 2: with performance (and terminal constraint)

            sysb_2 = obj.sys_block(diss.plant, X_f, G, diss.spec.rho);


            if diss.spec.target
                

                %optimize over the gain

                %hardcode the supply function

                rho = diss.spec.rho;
                rrecip = rho/(1-rho);
                gam = vars_spec.gam_p2p;
                M_u = -eye(diss.spec.nwp) * rrecip * (gam - mu);                

                               
                M_p2p = blkdiag(diss.iqc_rob.M, M_u);
                suppb_2 = obj.supply_block(plant_no_p, M_p2p);

                nzp = ssize(CDp, 1);

                %reciprocal by Schur complement
                M_yr =  (gam)*eye(nzp) *  (rrecip)^(-1);                

                con_M_2_corner = suppb_2 + sysb_2;

                %expand out the Schur complement
                con_M_2 = [con_M_2_corner, CDp'; CDp, M_yr];

                objective = gam;
            else
                
            
                
    
                %remember to take proper block diagonals and indexes!
                M_p2p = obj.merge_spec_M(diss.iqc_rob, diss.spec, vars_spec);
                
                
                suppb_2 = obj.supply_block(diss.plant, M_p2p);
                
                con_M_2 = sysb_2 + suppb_2;

                objective = 0;
            end

            %wrap it all together           

            sM1 = ssize(con_M_1,1);  sM2 = ssize(con_M_2,1);

            cons = append_lmi(cons, con_M_1 - obj.config.tol.M*eye(sM1), obj.LMILAB);   
            cons = append_lmi(cons, con_M_2 - obj.config.tol.M*eye(sM2), obj.LMILAB);                         

            con_M = {con_M_1, con_M_2};
        end


        function sol = process_recovery(obj, sol, lmi_out, alg_psi)
            %recover the controller

            %this code is with a full-order controller: duplication of the
            %number of internal model states
            %TODO: reduced order controller synthesis and recovery

            %get the system with the internal model
            P_trans =  obj.reg.connect_model(alg_psi, sol.rho);

            sys_cl = obj.system_closed_loop(P_trans, sol.vars.diss, sol.vars.reg, sol.vars.K);
            [Acl, Bcl, Ccl, Dcl] = ssdata(sys_cl);

            %for debugging
            G = obj.get_storage(sol.vars.diss, sol.vars.reg);

            %this is the (nonlinearly-warped) system that is certified as
            %possessing the desired performance and robustness
            %specifications
            sys_cal = ss(G \ Acl, G \ Bcl, Ccl, Dcl, 1);


            %evaluate the variables

            [A, B, C, D] = ssdata(P_trans);

            iz = [P_trans.index_z(), P_trans.index_zp()];
            iw = [P_trans.index_w(), P_trans.index_wp()];
            iu = P_trans.index_u();
            iy = P_trans.index_y();           

            nz = length(iz);
            nw = length(iw);
            nu = length(iu);
            ny = length(iy);

            Ak = sol.vars.K.A;
            Bk = sol.vars.K.B;
            Ck = sol.vars.K.C;
            Dk = sol.vars.K.D;
            
            K_warp = full([Ak, Bk; Ck, Dk]);

            [n] = size(Ak,1);
            % m = size(Bk);


            % S = (sol.vars.S);
            S = eye(n);

            Y = sol.vars.diss.GY;
            X = sol.vars.diss.GX;


            J = S - X * Y;
            [Up, Sig, Vp] = svd(J);

            % U = Up*Sig;
            ssig = sqrt(Sig);
            srsig = diag(1./(diag(ssig)));


            V = Vp*ssig;
            U = Up*ssig;

            Uinv = srsig*Up';
            Vinv = srsig*Vp';


            %similarity transformation


            %get right-side entries
            I = eye(n);
            Z1 = (Vinv*(I - X * Y')')';
            Z2 = (Vinv* (-U * Y')')';
            
            Z34 = [X, Z1; U, Z2] \ [zeros(n); eye(n)];
            
            Z3 = Z34(1:n, :);
            Z4 = Z34((n+1):end, :);
            
            T = [eye(n), Y'; zeros(n), V'];
            Ti = [eye(n), -Y' * Vinv'; zeros(n), Vinv'];
            
            SimG = [X, Z1; U, Z2];
            SimGi = [Y', Z3; V', Z4];

            %controller recovery

            Lblock = [Uinv, -Uinv*X*B(:, iu);
                zeros(nu, size(V, 2)), eye(nu)];

            LblockI = [U, X*B(:, iu); 
                zeros(nu, size(V, 2)), eye(nu)];

            Cblock = [Ak - X*A*Y, Bk;
                Ck, Dk];

            RblockI = [V' , zeros(size(V, 2), ny);
                C(iy, :)*Y, eye(ny)];

            Rblock = [Vinv', zeros(size(V, 2), ny);
                -C(iy, :)*Y*Vinv', eye(ny)];


            % Kblock0 = inv(LblockI)* (Cblock) * inv(RblockI);
            % Kblock1 = LblockI) \ Cblock) * inv(RblockI);
            % Kblockinv = RblockI') \ LblockI) \ Cblock)')';

            Kblock = Lblock * Cblock * Rblock;

            %extraction and exponential weighting
            Ac = Kblock(1:n, 1:n);
            Bc = Kblock(1:n, n+1:end);
            Cc = Kblock(n+1:end, 1:n);
            Dc = Kblock(n+1:end, n+1:end);

            K_nofeed_full = ss(Ac, Bc, Cc, Dc, 1);
            K_nofeed =minreal(K_nofeed_full,1e-5);

            


            %add the proper term by LFT
            D22 = D(iy, iu);
            Dfeed = zeros(nz+ny, nw+nu);            
            Dfeed(nz+1:end, nw+1:end) = D22;

            

            P_trans_nofeed = P_trans.ss;
            P_trans_nofeed.D = P_trans_nofeed.D - Dfeed;

            T_feed = [zeros(nu, ny), eye(nu); eye(ny), -D22];

            K_feed = lft(T_feed, K_nofeed, nu, ny);
            K_feed_full = lft(T_feed, K_nofeed_full, nu, ny);

            alg_trans = lft(P_trans, K_feed);
            alg_trans_nofeed = lft(P_trans_nofeed, K_nofeed);

            K_sub= rhotrafo(K_feed, 1/sol.rho);
            K_sub_full = rhotrafo(K_feed_full, 1/sol.rho);
            %connect the internal model: form the controller

            model = obj.reg.get_model(sol.vars.reg);

            K = lft(model, K_sub);
            K_full = lft(model, K_sub_full);

            
            %form the algorithm
            alg = lft(obj.sys.P, K);
            alg_full = lft(obj.sys.P, K_full);


            sol.alg = alg_full;
            sol.K = K_full;
            sol.model = model;
            sol.K_sub = K_sub;

            sol.alg_trans = alg_trans;  

            sol.gain = obj.validate_recovery_gain(alg_trans, sol.iqc_op_all);


            %verify performance of the algorithm
            %TODO: a postprocessing LMI (?) to check that the recovered 
            %algorithm satisfies the performance specifications 
        end

        function gain = validate_recovery_gain(obj, alg_trans, iqc_op_all)
            %VALIDATE_RECOVERY validate that the system obeys the stability
            %constraint (TODO: performance specs)


            %closed-loop and weighted system
            P = alg_trans.P(obj.sys.P.index_z, obj.sys.P.index_w);


            M = iqc_op_all.iqc.M;
            M = (M + M')/2;
            nw = floor(size(M, 1)/2);

            M11 = M(1:nw, 1:nw);
            M12 = M(nw + (1:nw), 1:nw);
            M22 = M(nw + (1:nw), nw + (1:nw));
            %is the constraint passive?
            is_passive = (norm(M11) + norm(M22) + norm(M12 - eye(nw)))==0;
            is_hinf = (norm(M11-eye(nw)) + norm(M22+eye(nw)) + norm(M12))==0;


            if is_passive
                gain_passive = -getPassiveIndex(-P, 'input');

                E=eye(nw);
                Tinf=[E sqrt(2)*E;sqrt(2)*E E];
                P_inf = lft(Tinf,P,nw,nw);

                gain_inf = norm(P_inf, 'inf');
            elseif is_hinf
                gain_inf = norm(P, 'inf');

                E=eye(nw);
                Tpass = [-E sqrt(2)*E;sqrt(2)*E -E];
                Ppass = lft(Tpass,P,nw,nw);

                gain_passive = -getPassiveIndex(-Ppass, 'input');
                elsep
                %TODO: advanced validation
                error('Customized validation is not yet implemented')
                gain_inf = 0;
                gain_passive = 0;
            end

            gain = [gain_passive, gain_inf];            

        end

    end

    
end


