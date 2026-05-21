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
            n = ssize(GX, 1);
            vars_diss= struct('GX', GX, 'GY', GY, 'GS', eye(n));

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

            [sys_cl, U_cl, V_cl] = obj.system_closed_loop(P, vars.diss, vars.reg, vars.K);
            
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
            [dyn_b, U_outer, V_outer] = obj.dynamics_block(sys_cl, quad);
            
            %wrap it all together
            objective = 0;

            con_M = stor_b + supp_b + dyn_b;


            if obj.elimination
                %knock them out

                V_elim = V_cl * V_outer;
                U_elim = U_cl * U_outer;
                

                U_null = null(U_elim, 'rational');
                V_null = null(V_elim, 'rational');

                con_M_U = U_null' * con_M * U_null;
                con_M_V = V_null' * con_M * V_null;

                sMU = ssize(con_M_U,1);
                sMV = ssize(con_M_V,1);

                cons = append_lmi(cons, con_M_U - obj.config.tol.M*eye(sMU), obj.LMILAB); 
                cons = append_lmi(cons, con_M_V - obj.config.tol.M*eye(sMV), obj.LMILAB); 

                %store the data
                con_M_0 = con_M;

                con_M = struct;
                %the main attributes for matrix elimination
                con_M.M0 = con_M_0;
                con_M.U = U_elim;
                con_M.V = V_elim;

                %subsidiaries for error checking
                con_M.U_null = U_null;
                con_M.V_null = U_null;
                con_M.M_U = con_M_U;
                con_M.M_V = con_M_V;

            else
                sM = ssize(con_M,1);
                cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 
            end
            %impose sign constraint            
            cons = obj.con_terminal(G, cons, [], diss.iqc_rob);
        end        


        function vars_new = augment_vars(obj, vars, diss, con_M)
            %AUGMENT_VARS add new variables/terms for recovery (useful for 
            %matrix elimination)                      
            vars_new = vars;            
            if obj.elimination
                vars_new.elim = con_M;            
            end

        end

        function [Ak, Ck] = recover_Ak_Ck(obj, vars_rec)
            %recover the Ak and Ck matrices
            %overridden by matrix elimination
            if obj.elimination
                nxi = size(vars_rec.K.B, 1);

                M0 = vars_rec.elim.M0;
                U = vars_rec.elim.U;
                V = vars_rec.elim.V;

                Un = null(U, 'rational');
                Vn = null(V, 'rational');


                con_U = Un' * M0 * Un;
                con_V = Vn' * M0 * Vn;
                AC_block = basiclmi(-M0, -U, V, 'Xmin');
                % AC_block = basiclmi(-M0, -U, V);

                Ak = AC_block(1:nxi, :);
                Ck = AC_block((nxi + 1):end, :);
            else
                Ak = vars_rec.K.A;
                Ck = vars_rec.K.C;
            end
        end

        function el = elimination(obj)
            %ELIMINATION is the matrix elimination lemma used?
            el = obj.config.syn.elimination;               
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


        function [sys_cl, U_cl, V_cl] = system_closed_loop(obj, P,  vars_diss, vars_reg, vars_K);
            %SYSTEM_CLOSED_LOOP closed-loop matrix after nonlinear
            %transformation

            %allow for matrix elimination
            %elimination: get rid of the [Ak; Ck] variables. 
            %solve only over [Bk; Dk].

            if obj.elimination
                %knock out the terms

                %get the variables
                GX = vars_diss.GX;
                GY = vars_diss.GY;

                
                Bk = vars_K.B;                
                Dk = vars_K.D;
                nk = ssize(Bk, 1);

                %should be a genplant type
                % P_net = diss.plant;


                [A, B, C, D] = ssdata(P);

                
                iu = P.index_u;
                iw = [P.index_w, P.index_wp];
               
                iy = P.index_y;
                iz = [P.index_z, P.index_zp];

                nxn = size(A, 1);
                nxi = ssize(Bk, 1);
                nu = length(iu);
                ny = length(iy);
                nz = length(iz);
                nw = length(iw);


                


                %closed loop without [Ak; Ck]
                Acal = [A*GY,  A + B(:, iu)*Dk*C(iy, :);
                    zeros(nk), GX*A + Bk*C(iy, :)];
                Bcal = [B(:, iw) + B(:, iu)*Dk*D(iy, iw);
                    GX*B(:, iw) + Bk*D(iy, iw)];
                Ccal = [C(iz, :)*GY, C(iz, :) + D(iz, iu)*Dk*C(iy, :)];
                Dcal = D(iz, iw) + D(iz, iu)*Dk*D(iy, iw);                

                sys_cl = sdpss(Acal, Bcal, Ccal, Dcal);

                
                %outer factors
                U_cl = [zeros(nxn, nxi), B(:, iu);
                    eye(nxi), zeros(nxi, nu);
                    zeros(nz, nxi), D(iz, iu)]';

                V_cl = [eye(nxi), zeros(nxi, nxn+nw)];


                %sys_cl(A; C) = sys_cl + U_cl [Ak; Ck] V_cl';
            else
                U_cl = [];
                V_cl = [];
                sys_cl = system_closed_loop@lmi_synthesis_interface(obj, P,  vars_diss, vars_reg, vars_K);
            end

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
            
        end


        

    end

    
end


