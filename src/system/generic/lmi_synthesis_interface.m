classdef lmi_synthesis_interface < lmi_dispatch_interface
    %LMI_ANALYSIS_INTERFACE 
    %Linear Matrix Inequality constraints for analysis of algorithmic
    %interconnections.
    %
    %
    %this is overridden by specialized analysis routines for system types:
    %   lti
    %   periodic
    %   switched robust
    %   switched jump
    
    properties
        reg; %internal model of the controller
                
        
    end

    
    methods
        function obj = lmi_synthesis_interface(sys, config)
            %LMI_SYNTHESIS_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_dispatch_interface(sys, config);


            %form the internal model
            reg_name = ['regulator_', sys.get_type()];

            reg_handle = str2func(reg_name);
            obj.reg = reg_handle(sys);

            %TODO: better options handling down below
            obj.config = config;
        end

        
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
            %
            %Output:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            
            
            
            %need to look up the right constraint            

            %Upper-levels: iterate over the systems

            if obj.elimination && (length(diss) > 1)
                error('Matrix Elimination (opt_config.syn.elimination=true) cannot be used if there is more than performance specification');
            end

            [cons, objective, vars, con_M] = cons_dynamic@lmi_dispatch_interface(obj, vars, cons, diss);

            

            %add new variables/terms for recovery (useful for matrix
            %elimination)

            
            vars = obj.augment_vars(vars, diss, con_M);
            
                      
        end


        %% variable creation

        function [vars, cons] = create_vars(obj, vars, cons, alg_psi, specs)
            %CREATE_VARS create the variables for the problem

            [vars.reg]  = obj.create_vars_regulator();
            [vars.diss, cons] = obj.create_vars_storage(cons, alg_psi);
            [vars.spec, cons] = obj.create_vars_spec(cons, specs);            
            [vars.K, cons]    = obj.create_vars_controller(cons, alg_psi);
        end

        function vars_new = augment_vars(obj, vars, diss, con_M)
            %AUGMENT_VARS add new variables/terms for recovery (useful for 
            %matrix elimination)
            vars_new = vars;            
        end

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

        function [vars_reg] = create_vars_regulator(obj)
            %CREATE_VARS_REGULATOR
            %parameterize the solutions to the regulator equations
            %use this as a variable in reduced-order control
            %
            %systems with more outputs than oracles can have freedom in the            
            %regulator equations (such as optimization problems with known 
            % Laplacian matrices)
            
            vars_reg = obj.reg.create_vars();            

        end


        function [GX, GY, cons] = define_storage_G(obj, cons, alg_psi, name)
            %DEFINE_STORAGE_G storage function for a specific subsystem

            %without terminal cost:
            n = ssize(alg_psi.A, 1);
            ns = obj.reg.ns;

            nX = n + ns;
           
            [GX, cons] = define_storage_GX(obj, cons, alg_psi, name);
            [GY, cons] = define_storage_GY(obj, cons, alg_psi, name);
            
        end


        function [GX, cons] = define_storage_GX(obj, cons, alg_psi, name);
             %DEFINE_STORAGE_GX primal storage function for a specific subsystem
             n = ssize(alg_psi.A, 1);
             ns = obj.reg.ns;

             nX = n + ns;
             GX = lmim(['GX', name], nX, nX, 'sym');
          

             if obj.config.tol.G_max < Inf                   
                cons = append_lmi(cons, obj.config.tol.GX_max*eye(nX)  - GX, obj.config.LMILAB);                            
                cons = append_lmi(cons, -obj.config.tol.G_min*eye(nX)  + GX, obj.config.LMILAB);                            
            end 
        end

        function [GY, cons] = define_storage_GY(obj, cons, alg_psi, name);
            %DEFINE_STORAGE_GY dual storage function for a specific subsystem            
            n = ssize(alg_psi.A, 1);
            ns = obj.reg.ns;
            
            nY = obj.get_GY_dim(n, ns);

            if nY > 0
                GY = lmim(['GY', name], nY, nY, 'sym');
            

    
                if obj.config.tol.G_max < Inf                   
                    cons = append_lmi(cons, obj.config.tol.GY_max*eye(nY)  - GY, obj.config.LMILAB);                
                    cons = append_lmi(cons, -obj.config.tol.G_min*eye(nY)  + GY, obj.config.LMILAB);                
                end 
            else
                GY = [];
            end

        end

        function ys = get_GY_dim(obj, n, ns)
            %dimension of the GY term
            ys = n +ns;
        end

        function verdict = reduced_order(obj)
            verdict = obj.config.syn.reduced_order;
        end

        function P_model = connect_model(obj, diss)

            P_model = obj.reg.connect_model(diss.plant, diss.rho);
        end

        function G = get_storage(obj, vars_diss, vars_reg)
            %GET_STORAGE get the storage function matrix G

            % 
            %
            GX = vars_diss.GX;
            GY = vars_diss.GY;
            GS = vars_diss.GS;
            
            G = [GY, GS; GS', GX];                
            
        end


        function [vars_spec, cons] = create_vars_spec(obj, cons, specs)
            %CREATE_VARS_SPEC declare variables for the specifications

            %maybe put this somewhere else?
            %
            %right now the variables are in the (spec) object.
            nspec = length(specs);
            vars_spec = cell(nspec, 1);
            for i = 1:nspec
                [vars_spec{i}, cons] = specs{i}.create_vars(cons, [], obj.config);
            end           
        end


        function cons = con_spread_single(obj, cons, GX, GY)
            %CON_SPREAD_SINGLE increase numerical conditioning by separating the 
            %primal and dual blocks
            np = ssize(GX, 1);
            spr = obj.config.tol.spread+1;
            cons_PH = [GX, (spr)*eye(np); (spr)*eye(np), GY];
            cons = append_lmi(cons, cons_PH, obj.LMILAB);

        end


        function cons = con_spread(obj, cons, vars)
            %CON_SPREAD increase numerical conditioning by separating the 
            %primal and dual blocks
            %invoke this over multiple subsystems
            % if ~obj.config.syn.reduced_order
                cons = obj.con_spread_single(cons, vars.diss.GX, vars.diss.GY);
            % end
        end

        function el = elimination(obj)
            el = false;
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
            
            include_Dk1 = true;
            [ny, nu] = size(D_mask);
            % if obj.config.syn.reduced_order
            %     %TODO: not yet implemented
            %     nc = n;
            %     nC = ny;
            %     include_Dk1 = false;
            % else
                nC = ns + ny;
                nc = n + ns;
            % end

            % ny = obj.sys.P.ny;
            % nu = obj.sys.P.nu;

            

            
            %declare the variables
            vars_K = struct;
            %easy: ABC
            

            if obj.elimination
                vars_K.A = [];                

                if obj.config.syn.elimination_type == 2    
                    %remove all terms [Ak, Bk; Ck, Dk]
                    %using triangular elimination (in development)
                    %lemma 4 of https://arxiv.org/pdf/1305.1746
                    vars_K.B = [];  
                    vars_K.C = [];  
                    vars_K.D = [];  
                elseif obj.config.syn.elimination_type == 1                                    
                    %remove [Ak, Bk; Ck1, Dk1]    
                    vars_K.B = [];  
                    vars_K.C = lmim(['Ck', name], ny, nc);
                    include_Dk1 = false;
                    vars_K.D = obj.form_Dk(alg_psi, D_mask, [], include_Dk1);
                    kq = [vars_K.C, vars_K.D];
                else
                    %remove [Ak; Ck]
                    vars_K.C = [];
                    vars_K.B = lmim(['Bk', name], nc, ny);
                    vars_K.D = obj.form_Dk(alg_psi, D_mask, [], include_Dk1);
                    kq = [vars_K.B;            
                    vars_K.D];
                    cons= append_lmi(cons, obj.config.tol.K_max*eye(sum(kq.dim)) - [zeros(kq.dim(1)), kq; kq', zeros(kq.dim(2))], obj.LMILAB);

            
                end
            else
                
                vars_K.A = lmim(['Ak', name], nc, nc);
                vars_K.B = lmim(['Bk', name], nc, ny);    
                vars_K.C = lmim(['Ck', name], nC, nc);
                vars_K.D = obj.form_Dk(alg_psi, D_mask, [], include_Dk1);

                kq = [vars_K.A, vars_K.B;            
                    vars_K.C,  vars_K.D];
                cons= append_lmi(cons, obj.config.tol.K_max*eye(sum(kq.dim)) - [zeros(kq.dim(1)), kq; kq', zeros(kq.dim(2))], obj.LMILAB);

            end

            
            
            %TODO: better interface here: number of inputs
            
            %bound entries of the controllers
            

        end

        %% formation of the Dk matrix in controller synthesis
        function D_mask = get_D_mask(obj)
            %GET_D_MASK get the direct feedthrough terms

            %the sparsity-constrained term for internal model control            
            D_mask_0 = obj.config.syn.D_mask;


            if isempty(D_mask_0)
                D_mask_0 = tril(ones(length(obj.sys.bind)));
            end

            %WARNING: do a better conversion on the coordinate lifts
            c = obj.sys.op{1}.c;
            D_mask = kron(D_mask_0, ones(c));

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

            D_mask = obj.get_D_mask();

            ns = size(obj.reg.R, 2);
            nu = obj.sys.nu;
            ny = obj.sys.ny;
            
            K_mask = logical([ones(nu + nxi + ns, nxi), [D_mask; ones(nxi+ns, ny)]]);

        end

        function [U, V] = get_K_tri_basis(obj, nxi)
            %GET_K_TRI_BASIS get a basis for the Dk
            %elimination method. Break up the lower-triangular Dk factor
            %to apply Lemma 4 of https://arxiv.org/pdf/1305.1746

            
            USE_LAST = all(nxi~=0);
            %get the mask for the permuted controller matrix
            K_mask = obj.get_K_mask(nxi);
            sz = size(K_mask);
            
            %now break it down
            %find the lower triangular factors
            [sel, cs] = max( K_mask ==0, [], 2 );
            [uf, ff] = unique(cs, 'first');
            
            coords= [uf, ff];
            
            coord_first = coords(2:end, :);
            coord_last = coords(1, 2);
            
            h_first = diff([1; coord_first(:, 1)]);
            
            coord_shift = coords(:, 2);
            coord_shift(1) = 0;
            
            
            %store the identity indexers            
            nc = size(coord_first, 1);
            U = cell(nc+USE_LAST, 1);
            V = cell(nc+USE_LAST, 1);            
            for i = 1:nc  
                % U{i} = speye(sz(1) - coord_shift(i), sz(1));
                U{i} = [sparse(sz(1)- coord_shift(i), coord_shift(i)),  speye(sz(1) - coord_shift(i))];
                hi = h_first(i);
                V{i} = sparse(1:hi, sum(h_first(1:i-1)) + (1:hi), ones(hi, 1), hi, sz(2));    
            end
            
            % U{nc+1} = speye(sz(1) - coord_last+1, sz(1));

            if USE_LAST
                U{nc+1} = [sparse(sz(1)-coord_last+1, coord_last - 1),  speye(sz(1) - coord_last+1)];
                hi = sz(2) - sum(cellfun(@(n) size(n, 1), V(1:end-1))); %not ideal
                V{nc+1} = sparse(1:hi, sum(h_first) + (1:hi), ones(hi, 1), hi, sz(2));
            end
        end


        function [Dk] = form_Dk(obj, alg_psi, D_mask, name, include_Dk1)
            %FORM_Dk: lower triangular structure needed for the controller
            %need a better interface for the mask


            %also, maybe an object structure for the internal model?
            
            if nargin < 4
                name = [];
            end

            if nargin < 5
                include_Dk1 = true;
            end

            
            n = ssize(alg_psi.A, 1);

            s = length(obj.sys.bind);
            
            ny = obj.sys.P.ny;            
            ns = obj.reg.ns;


            c = obj.sys.op{1}.c;


           
            %the unconstrained term for the internal model control

            if include_Dk1
                Dk1_var = lmim(['Dk1', name], ns, size(D_mask, 2), 'full');
                Dk = Dk1_var;
            else
                Dk = [];
            end


            % D_mask = D_mask_0;
 
            nd2= nnz(D_mask);
            if nnz(D_mask) > 0
                Dk2_var = lmim(['Dk2', name], 1, nd2, 'full');
    
                %make sure that the Dc2 term of the subcontroller is
                %lower-triangular
                
                counter = 1;
                for i = 1:size(D_mask, 1)
                     if any(D_mask(i, :))
                        eind = find(D_mask(i, :));
                        ncc = length(eind);
                        Dvar_mat = sparse(counter + (1:ncc)-1, 1:ncc, ones(ncc, 1), nd2, ncc);

                        Dvar = Dk2_var * Dvar_mat;
                        
                        Din = sparse(1:ncc, eind, ones(ncc, 1), ncc, size(D_mask, 1));
                        Din_var = Dvar * Din;

                        Dk = [Dk; Din_var];
    
                        
                        % Dk2_curr = sparse(, 1:nc )
                        % vars.Dk(i+(nu-s), j) = Dk2_var * eind;
                        counter = counter+nnz(D_mask(i, :));
                    else
                        Dk = [Dk; zeros(1, size(D_mask, 2))];
                    end
                    % end
                    
                % end
                end
            else
                Dk = [Dk; zeros(size(D_mask))];
            end

        end

        %% terminal constraints        
        function [cons, con_X] = con_terminal(obj, G, cons,  alg_psi, iqc_op)
            %CON_TERMINAL
            %terminal cost constraint (nonnegativity for the storage function G)
            %coupled positivity if the IQC has a terminal cost

            %too many arguments taken here
            X = iqc_op.X;

            %TODO: allow for reduced-order control           

            %TODO: check that this is the right formula, specifically when
            %X is a non-PSD terminal cost

            %matrix dilation results

            

            nf = ssize(X);
            n = ssize(G, 1);

            Ef = [eye(nf); zeros(n-nf, nf)];


            X_f = Ef * X * Ef';
            con_X = G + X_f;

            sx = ssize(con_X, 1);
            cons = append_lmi(cons, con_X - eye(sx)*obj.config.tol.X, obj.LMILAB);

        end

        %% helper functions to construct LMIs

        function stor_b = storage_block(obj, sys_cl, quad, G_curr, G_next)
            %STORAGE_BLOCK form the storage block in a synthesis problem
            if nargin < 5
                G_next = G_curr;
            end
               
            
            [n, nw] = ssize(sys_cl.B);
            nz = ssize(sys_cl.D, 1);
            nt = ssize(quad.U, 1);

            outer_curr = [eye(n, n); zeros(nw, n); zeros(n, n); zeros(nt, n)];

            outer_next = [zeros(n, n); zeros(nw, n); eye(n); zeros(nt, n)];
            
            stor_b = -outer_curr * G_curr * outer_curr'; 
            stor_b = stor_b - outer_next* G_next * outer_next'; 
        end


        function [dyn_b_he, U_outer, V_outer] = dynamics_block(obj, sys_cl, quad)
            %DYNAMICS_BLOCK form the supply block in a quadratic objective
            % problem


            [n, nw] = ssize(sys_cl.B);
            nz = ssize(sys_cl.D, 1);
            nt = ssize(quad.U, 1);

            center_cl = [sys_cl.A, sys_cl.B;
                sys_cl.C, sys_cl.D];

            outer_cl_right= [[eye(n), zeros(n, nw);
                zeros(nw, n), eye(nw)], zeros(n+nw, n+nt)];

            outer_cl_left = [zeros(n), zeros(n, nz);
                zeros(nw, n), quad.S;
                eye(n), zeros(n, nz); %check the sign here
                zeros(nt, n), quad.T];



            % outer_cl_left = [zeros(n), zeros(n, nw);
            %     zeros(nw, n), quad.S';
            %     eye(n), zeros(n, nw);
            %     zeros(nt, n), quad.T'];
            % 
            % outer_cl_right= [[eye(n), zeros(n, nz);
            %     zeros(nz, n), eye(nz)], zeros(n+nz, n+nt)];




            dyn_b = outer_cl_left * center_cl * outer_cl_right; 
            dyn_b_he = dyn_b + dyn_b';

            U_outer = -outer_cl_left';
            V_outer = outer_cl_right;

        end

        function supp_b = supply_block(obj, sys_cl, quad)
            % SUPPLY_BLOCK form the supply block in a quadratic objective
            % problem

            
            [n, nw] = ssize(sys_cl.B);
            nz = ssize(sys_cl.D, 1);
            nt = ssize(quad.U, 1);
            
            outer_Q = [zeros(n, nw); eye(nw); zeros(n, nw); zeros(nt, nw)];

            
            supp_b = outer_Q * (quad.Q + obj.config.tol.input_diss*eye(ssize(quad.Q, 1))) * outer_Q';

            outer_U = [zeros(n, nt); zeros(nw, nt); zeros(n, nt); eye(nt, nt)];

            if nt
                supp_b = supp_b - outer_U * quad.U * outer_U';
            end

        end

        function sys_cl = system_closed_loop(obj, P,  vars_diss, vars_reg, vars_K);
            %SYSTEM_CLOSED_LOOP closed-loop matrix after nonlinear
            %transformation

            GX = vars_diss.GX;
            GY = vars_diss.GY;

            %should be a genplant type
            % P_net = diss.plant;


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
            Acal = [A*GY + B(:, iu)*Ck,  A + B(:, iu)*Dk*C(iy, :);
                Ak, GX*A + Bk*C(iy, :)];
            Bcal = [B(:, iw) + B(:, iu)*Dk*D(iy, iw);
                GX*B(:, iw) + Bk*D(iy, iw)];
            Ccal = [C(iz, :)*GY+ D(iz, iu)*Ck, C(iz, :) + D(iz, iu)*Dk*C(iy, :)];
            Dcal = D(iz, iw) + D(iz, iu)*Dk*D(iy, iw);


            sys_cl = sdpss(Acal, Bcal, Ccal, Dcal);
        end


        function [cons, objective, con_M] = e2e_target(obj, vars, cons, diss)
            %E2E_target: energy to energy gain

            %is a special case of quadratic performance
            [cons, objective, con_M] = obj.quad(vars, cons, diss);

            %maybe keep this separate for analysis?
            
        end

        %% common specification calls

        
        %function [cons, objective, con_M] = quad(obj, vars, cons, diss)
        %Quadratic performance (defined on a per-system basis)

        %% Controller Recovery
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
            dissend = struct('plant', alg_psi, 'rho', sol.rho);
            P_trans =  obj.connect_model(dissend);
            
            %evaluate the variables
            [sol] = obj.recover_subcontroller(alg_psi, P_trans, sol);
                      
            %verify performance of the algorithm
            %TODO: a postprocessing LMI (?) to check that the recovered 
            %algorithm satisfies the performance specifications 
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

            [K_nofeed] = recover_subcontroller_warp(obj, P_trans, vars_rec);

            model = obj.reg.get_model(vars_rec.reg);

            K_report = obj.K_alg_report(P_trans, K_nofeed, model, rho);
            
            sol.alg_trans = K_report.alg_trans;
            sol.alg = lft(obj.sys.P, K_report.K);
            sol.model = K_report.model;           
            sol.K= K_report.K;
            sol.K_sub = K_report.K_sub;

            sol.gain = obj.validate_recovery_gain(sol.alg_trans, sol.iqc_op_all);

            
        end

        function [Ak, Bk, Ck, Dk] = recover_K_from_elim(obj, vars_rec)
            %recover the Ak and Ck matrices
            %overridden by matrix elimination
            Ak = vars_rec.K.A;
            Bk = vars_rec.K.B;            
            Ck = vars_rec.K.C;
            Dk = vars_rec.K.D;
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


            [A, B, C, D] = ssdata(P_trans);

            iz = [P_trans.index_z(), P_trans.index_zp()];
            iw = [P_trans.index_w(), P_trans.index_wp()];
            iu = P_trans.index_u();
            iy = P_trans.index_y();           

            nz = length(iz);
            nw = length(iw);
            nu = length(iu);
            ny = length(iy);
                        

            [Ak, Bk, Ck, Dk] = obj.recover_K_from_elim(vars_rec);                        

            S = (vars_rec.diss.GS);
            n = ssize(Ak, 1);

            Y = vars_rec.diss.GY;
            X = vars_rec.diss.GX;


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
        end


        function K_report = K_alg_report(obj, P_trans, K_nofeed, model, rho)
            %K_ALG_REPORT recover the algorithmic interconnection and the
            %controller
            
            D = P_trans.D;

            iz = [P_trans.index_z(), P_trans.index_zp()];
            iw = [P_trans.index_w(), P_trans.index_wp()];
            iu = P_trans.index_u();
            iy = P_trans.index_y();           

            nz = length(iz);
            nw = length(iw);
            nu = length(iu);
            ny = length(iy);

            %add the proper term by LFT
            D22 = D(iy, iu);
            Dfeed = zeros(nz+ny, nw+nu);            
            Dfeed(nz+1:end, nw+1:end) = D22;



            P_trans_nofeed = P_trans.ss;
            P_trans_nofeed.D = P_trans_nofeed.D - Dfeed;

            T_feed = [zeros(nu, ny), eye(nu); eye(ny), -D22];

            K_feed = lft(T_feed, K_nofeed, nu, ny);
            % K_feed_full = lft(T_feed, K_nofeed_full, nu, ny);


            alg_trans = lft(P_trans, K_feed);
            alg_trans_nofeed = lft(P_trans_nofeed, K_nofeed);


            K_sub= rhotrafo(K_feed, 1/rho);
            % K_sub_full = rhotrafo(K_feed_full, 1/rho);
            %connect the internal model: form the controller

            

            K = lft(model, K_sub);
            % K_full = lft(model, K_sub_full);


            % alg_full = lft(obj.sys.P, K_full);


            K_report = struct;
            
            K_report.K = K;
            K_report.model = model;
            K_report.K_sub = K_sub;

            K_report.alg_trans = alg_trans;  

        end

        function gain = validate_recovery_gain(obj, alg_trans, iqc_op_all)
            %VALIDATE_RECOVERY validate that the system obeys the stability
            %constraint (TODO: performance specs)


            %closed-loop and weighted system
            P = alg_trans.P(alg_trans.index_z, alg_trans.index_w);


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
            else
                %TODO: advanced validation
                warning('Customized validation is not yet implemented')
                gain_inf = 0;
                gain_passive = 0;
            end

            gain = [gain_passive, gain_inf];            

        end


    end

    methods (Abstract)
        %variable creation routines        
        quad(obj, vars, cons, diss)               
    end
end

