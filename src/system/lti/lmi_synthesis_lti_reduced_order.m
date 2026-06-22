classdef lmi_synthesis_lti_reduced_order < lmi_synthesis_lti
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
    %       quad
    %       stability
    %       e2e
    %
    %   TODO:
    %       p2p
    %       h2      
    %       e2p
    %
    %
    % Also need to figure out recovery
    
    methods
        function obj = lmi_synthesis_lti_reduced_order(sys, config)
            %LMI_SYNTHESIS_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_synthesis_lti(sys, config);
        end       
        
        %% reduced-order control indexers
        function [vars, cons] = create_vars(obj, vars, cons, alg_psi, specs)
            %CREATE_VARS create the variables for the problem

            [vars, cons]  = create_vars@lmi_synthesis_interface(obj, vars, cons, alg_psi, specs);

           
            vars.diss.GS = obj.Pibar(vars.diss, vars.reg);
           
        end

        function sol = process_recovery(obj, sol, lmi_out, alg_psi, diss)
            %recover the controller

            if nargin < 5
                diss = [];
            end
            
            %override this with other system types

            %this code is with a full-order controller: duplication of the
            %number of internal model states
            %TODO: reduced order controller synthesis and recovery

            %get the system with the internal model
            % dissend = struct('plant', diss{1}.plant_reg, 'rho', sol.rho);
            P_trans =  obj.connect_model(diss{1});
            
            %evaluate the variables
            [sol] = obj.recover_subcontroller(alg_psi, P_trans.P, sol);
                      
            %verify performance of the algorithm
            %TODO: a postprocessing LMI (?) to check that the recovered 
            %algorithm satisfies the performance specifications 
        end

        function [sol] = recover_subcontroller(obj, alg_psi, P_aug, sol)
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

            [K_nofeed, Gcl, Ycl] = recover_subcontroller_warp(obj, P_aug, vars_rec);


            
            model = obj.reg.get_model(vars_rec.reg);
            modelrho = rhotrafo(model, sol.rho);
            P_trans = lft(alg_psi, modelrho);

            K_report = obj.K_alg_report(P_trans, K_nofeed, model, rho);
            
            sol.alg_trans = K_report.alg_trans;
            sol.alg = lft(obj.sys.P, K_report.K);
            sol.model = K_report.model;           
            sol.K= K_report.K;
            sol.K_sub = K_report.K_sub;
            sol.gain = obj.validate_recovery_gain(sol.alg_trans, sol.iqc_op_all);

            sol.Gcl = Gcl;
            sol.Ycl = Ycl;
            
        end

        function P_model = connect_model(obj, diss)
            %CONNECT_MODEL connect the plant to the internal model

            %but this is reduced-order, so we play a different game here.

            % P_model is the augmented system following the loop
            % transformation and filter processing
            %
            % w  -> z
            % wp -> zp
            % d  -> e
            % u  -> y
            

            sys_aug = obj.sys;
            sys_aug.P = diss.plant_reg;

            [Plant, ~, alg_loop_aug] = sys_aug.build_plant(diss.iqc_data, diss.rho);

            P_model = struct('P', Plant, 'rho', diss.rho);
            % P_model = obj.reg.connect_model(diss.plant, diss.rho);
        end

        

        function Pb = Pibar(obj, vars_diss, vars_reg, invPi)
            %similarity transformation for optimization over Pi
            %used in regulator (reduced-order)  

            if nargin < 4
                invPi = false;
            end
            
            nxn = ssize(vars_reg.Pi, 1);
            ns = ssize(vars_reg.Pi, 2);
            
            nf = ssize(vars_diss.GX, 1) - nxn - ns;

            if invPi
                Pisign = 1;
            else
                Pisign = -1;
            end


            Pb = [[eye(nf), zeros(nf, nxn+ns)];
                   [zeros(nxn, nf), eye(nxn), Pisign*vars_reg.Pi]];

        end

        function [vars_K, cons] = create_vars_controller(obj, cons, alg_psi, name, D_mask)
            %CREATE_VARS_CONTROLLER create the nonlinearly-transformed
            %controller matrices

            %get the dimensions


            if nargin < 4
                name = [];
            end

            if nargin < 5
                D_mask = obj.get_D_mask;
            end

            n = ssize(alg_psi.A, 1);
            ns = obj.reg.ns;
            
            [ny, nu] = size(D_mask);


            %TODO: not yet implemented
            nc = n;
            nC = ny;
            include_Dk1 = false;


            % ny = obj.sys.P.ny;
            % nu = obj.sys.P.nu;


            
            %declare the variables
            vars_K = struct;
            %easy: ABC
            

            if obj.elimination
                vars_K.A = [];                

                
                %     %remove all terms [Ak, Bk; Ck, Dk]
                %     %using triangular elimination (in development)
                %     %lemma 4 of https://arxiv.org/pdf/1305.1746
                    vars_K.B = [];  
                    vars_K.C = [];  
                    vars_K.D = [];  
                % elseif obj.config.syn.elimination_type == 1                                    
                %     %remove [Ak, Bk; Ck1, Dk1]    
                %     vars_K.B = [];  
                %     vars_K.C = lmim(['Ck', name], ny, nc);                    
                %     vars_K.D = obj.form_Dk(alg_psi, D_mask, [], include_Dk1);
                %     kq = [vars_K.C, vars_K.D];
                % else
                %     %remove [Ak; Ck]
                %     vars_K.C = [];
                %     vars_K.B = lmim(['Bk', name], nc+ns, ny);
                %     vars_K.D = obj.form_Dk(alg_psi, D_mask, [], include_Dk1);
                %     kq = [vars_K.B;            
                %     vars_K.D];
                %     cons= append_lmi(cons, obj.config.tol.K_max*eye(sum(kq.dim)) - [zeros(kq.dim(1)), kq; kq', zeros(kq.dim(2))], obj.LMILAB);
                % 
                % 
                % end
            else
                
                if nc == 0
                    %static network and no filter dynamics: static output
                    %feedback subcontroller (before the internal model)
                    vars_K.A = zeros(nc+ns, nc);                    
                    vars_K.C = zeros(nC, nc);
                else
                    vars_K.A = lmim(['Ak', name], nc+ns, nc, 'full');                      
                    vars_K.C = lmim(['Ck', name], nC, nc, 'full');
                end
                vars_K.B = lmim(['Bk', name], nc+ns, ny, 'full');    
                vars_K.D = obj.form_Dk(alg_psi, D_mask, [], include_Dk1);

                kq = [vars_K.A, vars_K.B;            
                    vars_K.C,  vars_K.D];
                cons= append_lmi(cons, obj.config.tol.K_max*eye(sum(kq.dim)) - [zeros(kq.dim(1)), kq; kq', zeros(kq.dim(2))], obj.LMILAB);

            end

            
            
            %TODO: better interface here: number of inputs
            
            %bound entries of the controllers
            

        end

        function Ph = Pihat(obj, vars_diss, vars_reg, invPi)
            %similarity transformation for reduced-order control
            %used in regulator (reduced-order)

            if nargin < 4
                invPi = false;
            end
            nxn = ssize(vars_reg.Pi, 1);
            ns = ssize(vars_reg.Pi, 2);
            
            nf = ssize(vars_diss.GX, 1) - nxn - ns;

            
            Pb = obj.Pibar(vars_diss, vars_reg, invPi);
            Ph = [Pb; [zeros(ns, nf + nxn), eye(ns)]];
        end



        function cons = con_spread_single(obj, cons, GX, GY)
            %CON_SPREAD_SINGLE increase numerical conditioning by separating the 
            %primal and dual blocks
            % np = ssize(GX, 1);
            % spr = obj.config.tol.spread+1;           
            % cons_PH = [GX, (spr)*eye(np); (spr)*eye(np), GY];
            % cons = append_lmi(cons, cons_PH, obj.LMILAB);

            % cons = [];
        end

        function [sys_cl, U_cl, V_cl] = system_closed_loop(obj, Pr, vars_diss, vars_reg, vars_K);
            
            %acquire the transformed regulated expression
            
            GX = vars_diss.GX;
            GY = vars_diss.GY;

            
            P = Pr.P;

            rho = Pr.rho;

            [S, R] = obj.reg.exosystem();

            rhoi = (1/rho);

            Pibar = obj.Pibar(vars_diss, vars_reg);
            Pihat = obj.Pihat(vars_diss, vars_reg);
            Pihatinv = obj.Pihat(vars_diss, vars_reg, true);

            
            [A, B, C, D] = ssdata(P);
            iu = P.index_u;
            iw = [P.index_w];

            iy = P.index_y;
            iz = [P.index_z];


            ie = P.index_zp;
            id = P.index_wp;
            ns = size(S, 2);
            n = size(A, 1);

            % calligraphic matrices
            % from  convexification
            % [Y' Acl Y,  Y'Bcl ]
            % [Ccl Y,      Dcl  ]

            %follows formulation of (26) in 
            %https://www.sciencedirect.com/science/article/pii/S0005109808005402

            Ak = vars_K.A;
            Bk = vars_K.B;
            Ck = vars_K.C;
            Dk = vars_K.D;
            
            Aaug = [A, B(:, id); 
                zeros(ns, n), rhoi * S];

            Bpaug = [B(:, iw); 
                zeros(ns, length(iw))];
            Caug = [C(iy, :), D(iy, id)];

            ns = size(vars_reg.Pi, 2);
            nf = n - size(vars_reg.Pi, 1);
            
            Piaug = [zeros(nf, ns); vars_reg.Pi];
            Creg = [C(iw, :), -C(iw, :) * Piaug- D(iw, iu) * vars_reg.Gam];

            if obj.elimination

                if n ==0
                    %no filters nor network dynamics 
                    Acal = [zeros(size(Ak)), GX*Aaug];
                    Bcal = (GX*Pihatinv * Bpaug);
                    Ccal = (Creg);

                    
                else
                    %some filter or network dynamics
                    Acal = [A*GY,  Pibar * Aaug ;
                        zeros(n+nf+ns, n), GX*Aaug];
                    Bcal = [B(:, iw);
                        GX*Pihatinv * Bpaug];
                    Ccal = [C(iz, :)*GY, Creg];
                end
                Dcal = D(iz, iw);                

                %construct the outer factors

                nxn = size(B, 1);
                nxiU = n + ns;
                nu = length(iu);
                nz = length(iz);
                nw = length(iw);
                ny = length(iy);
                U_cl_base = [B(:, iu), zeros(nxn, nxiU);
                    zeros(nxiU, nu), eye(nxiU);
                    D(iz, iu), zeros(nz, nxiU)]';

                nxiV = n;
                V_cl_base = [eye(nxiV), zeros(nxiV, nxn + ns), zeros(nxiV, nw);
                    zeros(ny, nxiV), Caug, D(iy, iw)];

                %triangular decomposition
                [U_coord, V_coord] = obj.get_K_tri_basis([nxiU, nxiV]);

                ntri= length(V_coord);
                U_cl = cell(ntri+1, 1);
                V_cl = cell(ntri+1, 1);

                for i = 1:ntri
                    U_cl{i} = U_coord{i} * U_cl_base;
                    V_cl{i+1} = V_coord{i} * V_cl_base;
                end

                % error('reduced_order: elimination not yet supported')

                sys_cl = struct;
                sys_cl.A = Acal;
                sys_cl.B = Bcal;
                sys_cl.C = Ccal;
                sys_cl.D = Dcal;
            else
                if n ==0
                    %no filters nor network dynamics 
                    Acal = [Ak, GX*Aaug + Bk*Caug];
                    Bcal = (GX * Pihatinv * Bpaug + Bk*D(iy, iw));
                    Ccal = (Creg + D(iz, iu)*Dk*Caug);
                else
                    %some filter or network dynamics
                    Acal = [A*GY + B(:, iu)*Ck,  Pibar * Aaug + B(:, iu)*Dk*Caug;
                        Ak, GX*Aaug + Bk*Caug];
                    Bcal = [B(:, iw) + B(:, iu)*Dk*D(iy, iw);
                        GX * Pihatinv * Bpaug + Bk*D(iy, iw)];
                    Ccal = [C(iz, :)*GY+ D(iz, iu)*Ck, Creg + D(iz, iu)*Dk*Caug];
                end
                Dcal = D(iz, iw) + D(iz, iu)*Dk*D(iy, iw);

                U_cl = [];
                V_cl = [];
                sys_cl = sdpss(Acal, Bcal, Ccal, Dcal);
            end



            %this is a FORMAL system, `Acal' is rectangular and not square
            %due to the reduced-order structure
            



            % error('reduced_order: system closed loop block not yet supported')
        end

        function [Ak, Bk, Ck, Dk] = recover_K_from_elim(obj, vars_rec)
            %recover the eliminated matrices in the controller            
            if obj.elimination
                
                % https://www.sciencedirect.com/science/article/pii/0167691194000919
                %reconstruct the eliminated controller block
                M0 = vars_rec.elim.M0;
                U = vars_rec.elim.U;
                V = vars_rec.elim.V;                                

                
                    Knull = vars_rec.elim.null;

                    M_accum = M0;

                    ns = obj.reg.ns;
                    nu = obj.sys.nu;
                    ny = obj.sys.ny;

                    nxi1 = size(vars_rec.diss.GX, 1);
                    nxi2 = size(vars_rec.diss.GY, 1);

                    nxi = [nxi1; nxi2];
                    % K_block = zeros(ns + nxi(1) + nu, nxi(2) + ny);

                    [U_coord, V_coord] = obj.get_K_tri_basis(nxi);

                    K_block = zeros(size(U_coord{1}, 2), size(V_coord{1}, 2));
                    for i = (length(U)-1):-1:1
                        
                        %recover the current portion of the triangular
                        %controller
                        null_curr = Knull{i};
                        M_curr = null_curr' * M_accum * null_curr;
                        K_frag_curr = basiclmi(-M_curr, -U{i} * null_curr, V{i+1} * null_curr, 'Xmin');

                        K_embed_curr = U_coord{i}' * K_frag_curr * V_coord{i};
                        K_outer_curr = U{i}' * K_frag_curr * V{i+1};
                        K_block = K_block + K_embed_curr;

                        %prep for the next recovery step
                        M_accum = M_accum + K_outer_curr + K_outer_curr'; 
                    end

                    %now index the block

                    Ak = K_block(nu + (1:nxi(1)), 1:nxi(2));
                    Bk = K_block(nu + (1:nxi(1)), (nxi(2)+1):end);


                    Ck = K_block(1:nu, 1:nxi(2));
                    Dk = K_block(1:nu, (nxi(2)+1):end);

                    
            else
                [Ak, Bk, Ck, Dk] = recover_K_from_elim@lmi_synthesis_interface(obj, vars_rec);
            end
        end

        function K_mask = get_K_mask(obj, nxi)
            %K_mask: controller sparsity pattern
            %
            %[Ck2, Dk2
            % Ak,  Bk
            % Ck1, Dk1]
            %
            %Input: 
            %   nxi: number of controller states
            %Output:
            %   pattern K_mask

            %used for matrix elimination lemma for LTI systems

            if isscalar(nxi)
                nxi = nxi * [1, 1];
            end

            D_mask = obj.get_D_mask();

            ns = size(obj.reg.R, 2);
            nu = obj.sys.nu;
            ny = obj.sys.ny;
            
            K_mask = logical([ones(nu + nxi(1), nxi(2)), [D_mask; ones(nxi(1), ny)]]);

        end
        
        %% recovery

        function [K_nofeed, Xcal, Ycal] = recover_subcontroller_warp(obj, P_trans, vars_rec)

            %RECOVER_SUBCONTROLLER_WARP recover the nonlinearly warped
            %controller 
            %dynamics and indexers


            %for debugging
            % G = obj.get_storage(sol.vars.diss, sol.vars.reg);

            %this is the (nonlinearly-warped) system that is certified as
            %possessing the desired performance and robustness
            %specifications
            


            
            %get the regulator equation solution
            Pibar = obj.Pibar(vars_rec.diss, vars_rec.reg);
            Pihat= obj.Pihat(vars_rec.diss, vars_rec.reg);
            ihat = obj.Pihat(vars_rec.diss, vars_rec.reg, true);

            Gam = vars_rec.reg.Gam;
            Pi = vars_rec.reg.Pi;
            
            %get the base subcontroller

            %TODO: implement this
            [Ak, Bk, Ck, Dk] = obj.recover_K_from_elim(vars_rec);                        


            %get the closed-loop regulator equation
            Yred = vars_rec.diss.GY;
            X = vars_rec.diss.GX;

            Pibar = obj.Pibar(vars_rec.diss, vars_rec.reg);
            

            %transform back to the original coordinates
            % X = ihat' * Xhat * ihat;
            Xhat = X;

            n = ssize(Yred, 1);
            ns = obj.reg.ns;
            E = [eye(n); zeros(ns, n)]';
            Et = [zeros(ns, n), eye(ns)];

            if n==0
                Winv = [];
                Y = inv(X);
                W = [];
            else
                W = Yred - Pibar * inv(X) * Pibar';
                Winv = inv(W);
                Y = inv(X) + E' * W * E;
                
            end
            Z = Et*inv(X)*Pibar';


            %check the closed-loop system
            Ycal = [Yred, Pibar;
                   Z, [zeros(ns, n), eye(ns)];
                   -W', zeros(n, ns+n)];
            Ycalinv = inv(Ycal);
            



            %the similarity-transformed reduced-order storage matrix
            %what is solved in the program, should be positive definite (as
            %enforced by tolerances in config.tol)
            Gcl = obj.get_storage(vars_rec.diss, vars_rec.reg);

            %the original (full-order) storage matrix, should be positive
            %semidefinite (the reduced-order formalism relies on dropping
            %out modes)
            Gcl_large = [Y, eye(ns+n); eye(ns+n), X];



            Ut = X(:, 1:n);
            Vt = [-W; zeros(ns, n)];
            St = W;
            Ht = inv(St) + Ut'*inv(X)*Ut;

            ihat = inv(Pihat);


            %Xcal is the reduced-order storage matrices with inverse
            %Xcalinv. If they are not inverses of each other, then there is
            %a problem with the prior routines.            
            Xcal = [X, Ut;
                Ut', inv(St) + Ut'*(X \ Ut)];

            Xcalinv = [Y, Vt;
                Vt', inv(Ht) + Vt'*(Y \ Vt)];    



            Xhatcal = blkdiag(ihat', eye(n)) * Xcal * blkdiag(ihat, eye(n));
            Xhatcalinv = blkdiag(Pihat, eye(n)) * Xcalinv * blkdiag(Pihat', eye(n));  
            

            %verify the closed-loop behavior (for bugchecking) 
            cl_error = Ycal' * Xhatcal * Ycal - Gcl;



            diss_trans = struct('P', P_trans, 'rho', vars_rec.rho);
            sys_cl = obj.system_closed_loop(diss_trans, vars_rec.diss, vars_rec.reg, vars_rec.K);

            XAcl = Ycalinv' * sys_cl.A * Ycalinv;
            XBcl = Ycalinv' * sys_cl.B;
            Ccl = sys_cl.C * Ycalinv;
            Dcl = sys_cl.D;

            Acl = Xhatcalinv * XAcl;
            Bcl = Xhatcalinv * XBcl;



            %manual check of antipassivity
            nw = size(Dcl, 1);
            Ablock = [Acl, Bcl; eye(size(Acl)), zeros(2*n+ns, nw)];
            Sblock = kron([0, 1; 1, 0], eye(nw));
            Xblock = blkdiag(Xhatcal, -Xhatcal);
            Cblock = [Ccl, Dcl; zeros(nw, 2*n+ns), eye(nw)];
            ANTI = Ablock'*Xblock*Ablock + Cblock'*Sblock*Cblock;



            %this is the closed-loop response (after loop transformations
            % and multiplier augmentation), should satisfy the desired
            %performance specifications.
            sys_cl_rec = ss(Acl, Bcl, Ccl, Dcl, 1);

            A_rec = Acl(end-n+1 : end, end-n+1 : end); %this is the controller A matrix,             
            % should match with later recovery.

            %the performance of the recovered controller should match the
            %closed-loop quantity
            pass_rec = -getPassiveIndex(-sys_cl_rec, 'input');

            %now reconstruct a controller

            %recovery by transformation (preferred) or by solving a second
            %LMI with the given storage function Xcal (a bailout option)

            %index the system
            [A, B, C, D] = ssdata(P_trans);

            iu = P_trans.index_u;
            iw = [P_trans.index_w];
            iy = P_trans.index_y;
            iz = [P_trans.index_z];

            ie = P_trans.index_zp;
            id = P_trans.index_wp;        

            nz = length(iz);
            nw = length(iw);
            nu = length(iu);
            ny = length(iy);

            rhoi = 1/vars_rec.rho;

            %augmented system from the regulation conditions
            [S, R] = obj.reg.exosystem();
            Aaug = [A, B(:, id); zeros(ns, n), rhoi * S];            
            Bpaug = [B(:, iw); zeros(ns, length(iw))];
            Caug = [C(iy, :), D(iy, id)];
                        

            %original matrices in the system
            Agam = [A, -B(:, iu)*Gam; zeros(ns, n), rhoi * S];
            
            
            
            nf = n - size(Pi, 1);

            Piaug = [zeros(nf, ns); Pi];


            %perform controller recovery




            % 
            % %The transformation recovery is causing many issues.
            % %but we have the closed-loop response already. Use that:
            % %subspace arguments to extract the controller.
            % 
            % %having trouble finding such a realization. Maybe there's a
            % %coordinate transformation matrix occurring?
            % 
            % reg2 = obj.reg;
            % reg2.Gam = Gam;
            % reg2.Pi = Pi;
            % reg2.Phi = reg2.compute_Phi(Pi, Gam, []);
            % 
            % P_orc = P_trans.drop_performance;
            % 
            % P_model = reg2.connect_model(P_orc, 1/rhoi);
            % 
            % iu0 = P_model.index_u();
            % iu1 = iu0(1:ns);
            % iu2 = iu((ns+1):end);
            % 
            % nxn = n + ns;
            % nxi = n;
            % 
            % base_left = [zeros(nxn, nxi), P_model.Bu;
            %      eye(nxi), zeros(nxi, nu + ns);
            %     zeros(nz, nxi), P_model.Dzu];
            % 
            % base_right = [zeros(nxi, n + ns),  eye(nxi), zeros(nxi, nw);
            %     P_model.Cy, zeros(ny, nxi), P_model.Dyw];
            % 
            % Omega_base = [P_model.A, zeros(nxn, nxi), P_model.Bw;
            %     zeros(nxi, nxn), zeros(nxi), zeros(nxi, nw);
            %     P_model.Cz, zeros(nw, nxi), P_model.Dzw];
            % 
            % 
            % Omega_cl = [Acl, Bcl; Ccl, Dcl] - Omega_base;
            % 
            % Mat_sys = kron(base_right', base_left);
            % ans_sys = reshape(Omega_cl, [], 1);
            % 
            % 
            % %TODO: enforce the D_mask constraint in the system realization           
            % [ii, jj] = find(~obj.get_D_mask());
            % ND0 = length(ii);
            % 
            % ind_D2 = sub2ind([n+nu+ns, n+ny], n+ns+ii, n+jj);
            % Nmasked = length(ind_D2);
            % 
            % Mat_D = sparse(1:Nmasked, ind_D2, ones(Nmasked, 1), Nmasked, (n+nu+ns)*(n+ny));
            % ans_D = zeros(Nmasked, 1);
            % 
            % Mat_all = [Mat_sys; Mat_D];
            % ans_all = [ans_sys; ans_D];
            % 
            % %solve for a system in the subspace and verify
            % Kblock_vec = lsqminnorm(Mat_all, ans_all);
            % Kblock_vec(ind_D2) = 0; %make sure the zero entries are zero.
            % Kblock = reshape(Kblock_vec, n + length(iu0), n + ny);
            % 
            % 
            % rec_error = norm(Mat_all*Kblock_vec - ans_all);

            

            %(33-34) of https://www.sciencedirect.com/science/article/pii/S0005109808005402

            if n==0
                Akt = X \ Ak;                
            else
                Akt = X \ Ak - Aaug * Y * Pibar';

            end
            Bkt = X \ Bk;
            Ckt = Ck + Gam*Z;
            Dkt = Dk;

            LblockI = [eye(n), rhoi * Pi, B(:, iu);
                zeros(ns, n), rhoi * eye(ns), zeros(ns, nu);
                zeros(nu, n), zeros(nu, ns), eye(nu)];

            Cblock = [Akt, Bkt;
                Ckt, Dkt];


            % RblockI = [-W, zeros(size(Winv, 1), ny);
            %       Caug * Y * Pibar', eye(ny)];


            %recovery conditions involve Caug, not Creg.
            Rblock = [-Winv, zeros(size(Winv, 1), ny);
                Caug * Y * Pibar' * Winv, eye(ny)];

            Kblock = LblockI \ (Cblock * Rblock);

            %extraction and exponential weighting


            Ac = Kblock(1:n, 1:n);
            Bc = Kblock(1:n, n+1:end);
            Cc = Kblock(n+1:end, 1:n);
            Dc = Kblock(n+1:end, n+1:end);



            K_nofeed_full = ss(Ac, Bc, Cc, Dc, 1);
            K_nofeed = K_nofeed_full;            

            % if norm(rec_error) > 1e-6
            %     error('Reduced Order LTI: failure of subspace-based controller reconstruction')
            % end
        end



        

    end

    
end


