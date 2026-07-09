classdef lmi_synthesis_periodic_orbit < lmi_synthesis_interface
    %LMI_SYNTHESIS_PERIODIC_ORBIT synthesisLMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
    %
    % Orbit structure on the periodicity
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


    methods
        function obj = lmi_synthesis_periodic_orbit(sys,config)
            %LMI_SYNTHESIS_PERIODIC undefined
            %   undefined
            obj@lmi_synthesis_interface(sys, config);
        end

        %% definition of variables and helpers

        function ns = Nss(obj)
            %NSS: Number of subsystems            
            ns = obj.sys.Nss;
        end

        %% main call

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
            n = ssize(GX, 1);

            GS = eye(n);
            vars_diss= struct('GX', GX, 'GY', GY, 'GS', GS);

        end

        function P_model = connect_model(obj, diss)

            P_model = obj.reg.connect_model(diss.plant, 1, diss.rho);
        end

        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance

            %TODO: figure out why this doesn't work.


            %get the variables of the problem            
            % Gnext = obj.get_storage(vars.diss, vars.reg);
            Gcurr = obj.get_storage(vars.diss, vars.reg);

            n = ssize(vars.diss.GX, 1);
            R = obj.sys.R;
            c = size(R, 1);

            Rkron = kron(eye(2*n/c), R);
            Rkroninv = inv(Rkron);

            Rkron_small = Rkron(1:n, 1:n);
            
            Gnext = Rkron' * Gcurr * Rkron;
            % Gnext = Gcurr;




            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            P = obj.connect_model(diss);

            %account for the periodicity
            vars_diss = vars.diss;
            vars_diss.GX = Rkron_small' * vars.diss.GX * Rkron_small;            

            [sys_cl, U_cl, V_cl] = obj.system_closed_loop(P, vars_diss, vars.reg, vars.K);
            
            %index the quadratic specification
            if isempty(diss.spec.izp)
                ind_p = 1:(diss.iqc_rob.nz);
                ind_q = diss.iqc_rob.nz + (1:(diss.iqc_rob.nw));
            else
                ind_p = 1:(diss.iqc_rob.nz + diss.spec.izp);
                ind_q = (diss.iqc_rob.nz + diss.spec.izp) + (1:(diss.iqc_rob.nw + diss.spec.iwp));
            end
            
            vars_spec = vars.spec{diss.spec.id};
            np = diss.iqc_rob.np;
            nq = diss.iqc_rob.nq;
            M_quad_rob = quad_objective_decomp(diss.iqc_rob.M, 1:np, np + (1:nq));

            [M_quad_spec, objective] = diss.spec.supply_quad(vars_spec);

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
            % con_M = 2;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 


            %impose sign constraint            
            cons = obj.con_terminal(Gcurr, cons, [], diss.iqc_rob);
        end

        %TODO: e2e_target

        %% Recovery
        % Call the recovery method to finalize the synthesis process
        % [cons, objective, con_M] = obj.recovery(vars, cons, diss);

        function gain = validate_recovery_gain(obj, alg_trans, iqc_op_all)
            %VALIDATE_RECOVERY validate that the system obeys the stability
            %constraint (TODO: performance specs)

            %use the monodromy system to get specs
            n = alg_trans.dump_dim();

            sys_trans = obj.sys;
            sys_trans.K = [];
            sys_trans.P = alg_trans;
            sys_trans.P.P;
            %TODO: fix this
            alg_trans_lti = periodic_lift(sys_trans);


            gain = validate_recovery_gain@lmi_synthesis_interface(obj, alg_trans_lti, iqc_op_all);
        end

        
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
            [K_nofeed, Gcl, Ycl] = recover_subcontroller_warp(obj, P_trans, vars_rec);


            %package it up
            
            model = obj.reg.get_model(1, vars_rec.reg);

            K_report = obj.K_alg_report(P_trans, K_nofeed, model, rho);

            %form the algorithm
            sol.cert.alg_trans = K_report.alg_trans;
            sol.cert.alg = lft(obj.sys.P, K_report.K);
            sol.cert.model = K_report.model;           
            sol.K= K_report.K;
            sol.cert.K_sub = K_report.K_sub;                

        
            sol.cert.Gcl = Gcl;
            sol.cert.Ycl = Ycl;


            sol.gain = obj.validate_recovery_gain(sol.cert.alg_trans, sol.cert.iqc_op_all);
        end


        function [K_nofeed, Gcal, Ycal] = recover_subcontroller_warp(obj, P_trans, vars_rec)

            %RECOVER_SUBCONTROLLER_WARP recover the nonlinearly warped
            %controller 
            %dynamics and indexers


            %for debugging
            

            %this is the (nonlinearly-warped) system that is certified as
            %possessing the desired performance and robustness
            %specifications
            % sys_cl = obj.system_closed_loop(P_trans, sol.vars.diss, sol.vars.reg, sol.vars.K);

            % sys_cal = ss(G \ Acl, G \ Bcl, Ccl, Dcl, 1);


            %extract the storage variables and the factorizations


            K_nofeed = cell(obj.Nss, 1);        

                Y = vars_rec.diss.GY;
                X = vars_rec.diss.GX;

                S = vars_rec.diss.GS;

                J = S - X * Y;
                [Up, Sig, Vp] = svd(J);

                % U = Up*Sig;
                ssig = sqrt(Sig);
                srsig = diag(1./(diag(ssig)));


                V = Vp*ssig;
                U = Up*ssig;

                Uinv = srsig*Up';
                Vinv = srsig*Vp';

                n = size(X, 1);
                c = size(obj.sys.R, 1);
                Rkron = kron(eye(n/c), obj.sys.R);
                Rkroninv = inv(Rkron);
                Uinvnext = Rkroninv' * Uinv * Rkroninv; % Store the current Vinv for the next iteration
                Xnext = Rkron' * X * Rkron; % Update Yprev for the next iteration
            
            %get the indexers
            Pt = P_trans;            

            iz = [Pt.index_z(), Pt.index_zp()];
            iw = [Pt.index_w(), Pt.index_wp()];
            iu = Pt.index_u();
            iy = Pt.index_y();           

            nz = length(iz);
            nw = length(iw);
            nu = length(iu);
            ny = length(iy);

            %the storage matrices and transformation

            G = obj.get_storage(vars_rec.diss, vars_rec.reg);

            Ycal = [Y, eye(n); V', zeros(n)];
            iYcal = inv(Ycal);
            Gcal = iYcal' * G * iYcal;


            
                [A, B, C, D] = ssdata(P_trans);
                Ak = vars_rec.K.A;
                Bk = vars_rec.K.B;
                Ck = vars_rec.K.C;
                Dk = vars_rec.K.D;
    
                
                n = ssize(Ak, 1);
                
    
    
                %controller recovery
    
                Lblock = [Uinvnext, -Uinvnext*Xnext*B(:, iu);
                    zeros(nu, size(V, 2)), eye(nu)];
    
                
                Cblock = [Ak - Xnext*A*Y, Bk;
                    Ck, Dk];
  
    
                Rblock = [Vinv', zeros(size(Vinv, 2), ny);
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
    
                K_nofeed= ss(Ac, Bc, Cc, Dc, 1);                
            
        end
    end
end