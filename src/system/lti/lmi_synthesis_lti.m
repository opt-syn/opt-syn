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

        function sys_cl = system_closed_loop(obj, vars_diss, vars_reg, vars_K, diss);
            %SYSTEM_CLOSED_LOOP closed-loop matrix after nonlinear
            %transformation

            GX = vars_diss.GX;
            GY = vars_diss.GY;

            %should be a genplant type
            % P_net = diss.plant;

            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            P = obj.reg.connect_model(diss.plant);

            [A, B, C, D] = ssdata(P);
            iu = P.index_u;
            iw = [P.index_w, P.index_wp];

            iy = P.index_y;
            iz = [P.index_z, P.index_zp];

            % calligraphic matrices
            % from  convexification
            % [Y' Acl Y,  Y'Bcl ]
            % [Ccl Y,      Dcl  ]


            Ak = vars_K.A;
            Bk = vars_K.B;
            Ck = vars_K.C;
            Dk = vars_K.D;
            %
            Acal = [A*GX + B(:, iu)*Ck,  A + B(:, iu)*Dk*C(iy, :);
                    Ak, GY*A + Bk*C(iy, :)];
            Bcal = [B(:, iw) + B(:, iu)*Dk*D(iy, iw);
                GY*B(:, iw) + Bk*D(iy, iw)];
            Ccal = [C(iz, :)*GX + D(iz, iu)*Ck, C(iz, :) + D(iz, iu)*Dk*C(iy, :)];
            Dcal = D(iz, iw) + D(iz, iu)*Dk*D(iy, iw);
    

            sys_cl = sdpss(Acal, Bcal, Ccal, Dcal);
        end
        



        %% Quadratic performance (infinite horizon)


        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance

            

            %get the variables of the problem
            G = obj.get_storage(vars.diss, vars.reg);
            sys_cl = obj.system_closed_loop(vars.diss, vars.reg, vars.K, diss);

            
            %index the quadratic specification
            vars_spec = vars.spec{diss.spec.id};
            M_quad = obj.merge_spec_M(diss.iqc_rob, diss.spec, vars_spec);

            if isempty(diss.spec.izp)
                ind_p = 1:(diss.iqc_rob.nz);
                ind_q = diss.iqc_rob.nz + (1:(diss.iqc_rob.nw));
            else
                ind_p = 1:(diss.iqc_rob.nz + diss.spec.izp);
                ind_q = (diss.iqc_rob.nz + diss.spec.izp) + (1:(diss.iqc_rob.nw + diss.spec.iwp));
            end
            %use eigenvalue arguments here

            Qq = M_quad(ind_p, ind_p);
            Sq = M_quad(ind_p, ind_q);
            Rq = M_quad(ind_q, ind_q);


            [RqV, RqD] = eig(Rq);
            eRq = diag(RqD);
            ind_pos = find(abs(eRq) > 1e-12);

            Tq = RqV(:, ind_pos);
            Uq = diag(1./eRq(ind_pos));
            
            %formulation from ParDynSyn notes (parametric dynamic
            %synthesis)

            %acquire the dimensions
            n = ssize(sys_cl.A, 1);
            nw = ssize(sys_cl.B, 2);
            nz = ssize(sys_cl.C, 1);
            nt = length(ind_pos);

            %[nx, nz, nx, nt]

            %TODO: audit this, break up into other routines
            outer_Q = [zeros(n, nz); eye(nz); zeros(n, nz); zeros(nt, nz)];
           
            supp_b = 0;
            supp_b = -outer_Q * Qq * outer_Q';

            outer_U = [zeros(n, nt); zeros(nz, nt); zeros(n, nt); eye(nt, nt)];

            if nt
                supp_b = supp_b + outer_U * Uq * outer_U';
            end
            
            %supply block building
            outer_curr = [diss.spec.rho*eye(n, n); zeros(nz, n); zeros(n, n); eye(nt, n)];

            outer_next = [zeros(n, n); zeros(nz, n); eye(n); eye(nt, n)];

            G_curr = G;
            G_next = G;
            sys_b_G = outer_curr * G_curr * outer_curr'; 
            sys_b_G = sys_b_G + outer_next* G_next * outer_next'; 


            %now for the controller parameters
            %TODO: verify dimensions here
            outer_cl_left = [zeros(n), zeros(n, nw);
                zeros(nw, n), Sq;
                eye(n), zeros(n, nw);
                zeros(nt, n), Tq];

            outer_cl_right= [[eye(n), zeros(n, nz);
                zeros(nz, n), eye(nz)], zeros(n+nz, n+nt)];
                

            center_cl = [sys_cl.A, sys_cl.B;
                sys_cl.C, sys_cl.D];

            dyn_b = outer_cl_left * center_cl * outer_cl_right; 
            dyn_b_he = dyn_b + dyn_b';
            % sys_b = sys_b - dyn_b_he;

            %wrap it all together
            objective = 0;

            con_M = sys_b_G + supp_b + dyn_b_he;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint
            %change this up
            cons = obj.con_terminal(G, cons, [], diss.iqc_rob);
        end        



        function [cons, objective, con_M] = e2e_target(obj, vars, cons, diss)
            %E2E_TARGET: use a Schur complement to minimize the energy to
            %energy gain of the transfer function

            G = vars.diss.G;

            
           
            sysb = obj.sys_block(diss.plant, G, G, diss.spec.rho);

            %variable to optimize
            mu = vars.spec{diss.spec.id}.mu_l2;

            [plant_no_p, CDp] = obj.separate_performance_output(diss);

            %form the supply
            nwp = length(diss.spec.iwp);
            M_base = blkdiag(diss.iqc_rob.M, -mu * eye(nwp));
            
            objective = mu;            
            suppb = obj.supply_block(plant_no_p, M_base);


            %wrap it all together           
            con_M_corner = sysb + suppb;
            nzp = ssize(CDp, 1);
            con_M = [con_M_corner, CDp'; CDp, mu*eye(nzp)];


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.tol.M*eye(sM), obj.LMILAB);   
            
            
            %impose sign constraint
            cons = obj.con_terminal(G, cons, diss.iqc_rob);
        end

        %% Peak-to-Peak norm (at each finite horizon)

        function [cons, objective, con_M] = p2p(obj, vars, cons, diss)
            %p2p: certificate of peak to peak induced norm
            %
            % sup norm(zp, 2) / norm(wp, 2) <= objective

            % verification by Theorem 4 of https://www.sciencedirect.com/science/article/pii/S2405896323008194

            %storage matrix
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
            
            sysb_1 = obj.sys_block(diss.plant, G, G, diss.spec.rho);
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

            cons = append_lmi(cons, con_M_1 - obj.tol.M*eye(sM1), obj.LMILAB);   
            cons = append_lmi(cons, con_M_2 - obj.tol.M*eye(sM2), obj.LMILAB);                         

            con_M = {con_M_1, con_M_2};
        end


        function sol = process_recovery(obj, sol, lmi_out, alg_psi)
            %recover the controller

            %this code is with a full-order controller: duplication of the
            %number of internal model states
            %TODO: reduced order controller synthesis and recovery

            %get the system with the internal model
            P_trans =  obj.reg.connect_model(alg_psi);

            %evaluate the variables

            [A, B, C, D] = ssdata(P_trans);

            iz = [P_trans.index_z(); P_trans.index_zp()];
            iw = [P_trans.index_w(); P_trans.index_wp()];
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

            G = obj.get_storage(sol.vars.diss, sol.vars.reg);

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

            Ac = Kblock(1:n, 1:n);
            Bc = Kblock(1:n, n+1:end);
            Cc = Kblock(n+1:end, 1:n);
            Dc = Kblock(n+1:end, n+1:end);

            K_sub_full = ss(Ac, Bc, Cc, Dc, 1);
            K_sub=minreal(K_sub_full,1e-5);

            


            %add the proper term by LFT
            %this part may be incorrect


            D22 = D(iy, iu);
            Dfeed = zeros(nz+ny, nw+nu);            
            Dfeed(nz+1:end, nw+1:end) = D22;

            

            P_trans_nofeed = P_trans.ss;
            P_trans_nofeed.D = P_trans_nofeed.D - Dfeed;

            T_feed = [zeros(nu, ny), eye(nu); eye(ny), -D22];

            K_feed = lft(T_feed, K_sub, nu, ny);
            K_feed_full = lft(T_feed, K_sub_full, nu, ny);

            %connect the internal model: form the controller

            model = obj.reg.get_model(sol.vars.reg);

            K = lft(model, K_feed);
            K_full = lft(model, K_feed_full);

            
            %form the algorithm
            alg = lft(obj.sys.P, K);
            alg_full = lft(obj.sys.P, K_full);


            sol.alg = alg_full;
            sol.K = K_full;
            sol.model = model;
            sol.K_sub = K_sub;

            sol.alg_trans = lft(P_trans, K_feed);  


            %verify performance of the algorithm
            %TODO: a postprocessing LMI (?) to check that the recovered 
            %algorithm satisfies the performance specifications 

            



        end

    end

    
end


