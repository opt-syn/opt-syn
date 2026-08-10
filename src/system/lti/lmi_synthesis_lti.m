classdef lmi_synthesis_lti < lmi_synthesis_interface
    %LMI_SYNTHESIS_LTI synthesis LMIs for algorithmic interconnections
    %involving linear-time-invariant (LTI) networks and controllers
    


    
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
    %       stability
    %       e2e
    %       quad
    %
    %   TODO:
    %       p2p
    %       h2      
    %       e2p

    
    methods
        function obj = lmi_synthesis_lti(sys, config)
            %LMI_SYNTHESIS_LTI Constructor
            obj@lmi_synthesis_interface(sys, config);


            if obj.config.syn.elimination
                %check D_mask, only perform elimination if it is
                %block-lower triangular
                obj.config.syn.elimination = obj.check_lower_triangular();
            end
            

        end       
        
        function [is_tri] = check_lower_triangular(obj)
            %CHECK_LOWER_TRIANGULAR is the D matrix constrained to be
            %block-lower-triangular? This must be true to use matrix
            %elimination
            %
            %Returns:
            %   is_tri (bool): verdict on lower triangularity

            D_mask = obj.get_D_mask();           
            is_tri = true;
            i0 = find(D_mask == 0);
            [i, j] = ind2sub(size(D_mask), i0);

            for k = 1:length(i)
                D_right = D_mask(i(k), (j(k)+1):end);
                D_top = D_mask(1:(i(k)-1), j(k));

                if sum(D_right) || sum(D_top)
                    is_tri = false;
                    break
                end
            end

        end
        
        %% definition of variables and helpers
        function [vars_diss, cons]= create_vars_storage(obj, cons, alg_psi, name)
            %create_vars_storage create variables for the dissipation
            %constraints
            %
            %Args:
            %   cons:       accumulated constraints
            %   alg_psi:    the filtered algorithmic interconnection
            %   name:       a name for the variable
            %Returns:
            %   vars_diss:   variables of the problem in the dissipation constraints
            %   cons:   accumulated constraints
            
            if nargin < 4
                name = [];
            end

            
            [GX, GY, cons] = obj.define_storage_G(cons, alg_psi, name);
            n = ssize(GX, 1);

            GS = eye(n);
            vars_diss= struct('GX', GX, 'GY', GY, 'GS', GS);

        end
        

        %% Quadratic performance (infinite horizon)        
        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance
            %
            %Args:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss (diss_data):   structure describing the dissipation constraint
            %Returns:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            
            %   con_M:      PSD blocks for the dynamics constraint

            %get the variables of the problem
            G = obj.get_storage(vars.diss, vars.reg);
            
            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            rhou = obj.used_rho(diss);
            P = obj.connect_model(diss, rhou);            
            
            [sys_cl, U_cl, V_cl] = obj.system_closed_loop(P, vars.diss, vars.reg, vars.K, diss.rho);
            
            %index the quadratic specification
            vars_spec = vars.spec{diss.spec.id};

            np = diss.iqc_rob.np;
            nq = diss.iqc_rob.nq;

            M_quad_rob = quad_objective_decomp(diss.iqc_rob.M, 1:np, np + (1:nq));
            [M_quad_spec, objective] = diss.spec.supply_quad(vars_spec);

            quad = obj.merge_quad(M_quad_rob, M_quad_spec);

            %the quadratic objective
            supp_b = obj.supply_block(sys_cl, quad);

            %the storage
            stor_b = obj.storage_block(sys_cl, quad, G, G);

            %the dynamics
            [dyn_b, U_outer, V_outer] = obj.dynamics_block(sys_cl, quad);
            
            %wrap it all together           
            con_M = -(stor_b + supp_b + dyn_b);


            if obj.elimination
                %knock them out

                if ~iscell(U_cl)
                    U_cl0 = U_cl;
                    V_cl0 = V_cl;
                    U_cl = {U_cl0, []};
                    V_cl = {[], V_cl0};
                end
                
                ntri = length(U_cl);                    

                U_elim = cell(ntri, 1);
                V_elim = cell(ntri, 1);
                con_M_null = cell(ntri, 1);
                V_accum = [];
                null_accum = cell(ntri, 1);
                for i = 1:ntri

                    if ~isempty(U_cl{i})
                        U_elim{i} = U_cl{i} * U_outer;
                    end
                    if ~isempty(V_cl{i})
                        V_elim{i} = V_cl{i} * V_outer;
                    end
                    V_accum = [V_accum; V_cl{i}];
                    
                    %get the nullspace
                    if isempty(V_accum)
                        elim_curr = U_elim{i}; %used in elimination
                        null_store = eye(ssize(con_M, 1));       %stored for recovery
                    else
                        elim_curr = [V_accum * V_outer; U_elim{i}];
                        null_store = null(V_accum * V_outer, 'rational');
                    end
                    null_accum{i} = null_store;

                    null_curr = null(elim_curr, 'rational');


                    %form and enforce the constraint
                    con_M_curr = null_curr'* con_M * null_curr;

                    sU = ssize(con_M_curr, 1);
                    con_M_null{i} = con_M_curr;
                    
                    cons = append_lmi(cons, con_M_curr - obj.config.tol.M*eye(sU), obj.LMILAB); 

                end

                con_M_0 = con_M;

                con_M = struct;
                %the main attributes for matrix elimination
                con_M.M0 = con_M_0;
                con_M.U = U_elim;
                con_M.V = V_elim;
                con_M.null = null_accum;
                con_M.Mnull = con_M_null;
            else
                sM = ssize(con_M,1);
                cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 
            end
            %impose sign constraint            
            cons = obj.con_terminal(G, cons, [], diss.iqc_rob);
        end        



        function ys = get_GY_dim(obj, n, ns)
            %dimension of the GY term
            %Args:
            %   n:  number of states
            %   ns: number of exogenous disturbances
            %Returns:
            %   ys: size of GY matrix
            if obj.reduced_order
                ys = n;
            else
                ys = n +ns;
            end
        end



       function G = get_storage(obj, vars_diss, vars_reg)
            %GET_STORAGE get the storage function matrix G
            %
            %Args:                   
            %   vars_diss:   variables of the problem in the dissipation constraints
            %   vars_reg:   variables for regulator equation
            %Returns:            
            %   G:   the closed-loop storage matrix (warped)
            
            GX = vars_diss.GX;
            GY = vars_diss.GY;
            
            GS = vars_diss.GS;
            
            G = [GY, GS; GS', GX];                
            
       end

       %% Peak-to-Peak norm (at each finite horizon)

       function [cons, objective, con_M] = p2p(obj, vars, cons, diss)
            %p2p: certificate of finite-horizon peak-to-peak norm bounds
            % when starting at a zero (steady state) initial condition, not
            % transient performance.
            %
            %Args:
            %   cons:       accumulated constraints
            %   alg_psi:    the filtered algorithmic interconnection
            %   name:       a name for the variable
            %Returns:
            %   vars_diss:   variables of the problem in the dissipation constraints
            %   cons:   accumulated constraints
            %
            %Warning:
            %   not yet stable. do not use yet.
            
            %get the variables of the problem
            G = obj.get_storage(vars.diss, vars.reg);
            
            %IMPORTANT!                 
            rhou = obj.used_rho(diss);
            P = obj.connect_model(diss, rhou);

            [sys_cl, U_cl, V_cl] = obj.system_closed_loop(P, vars.diss, vars.reg, vars.K);
            
            %index the quadratic specification
            vars_spec = vars.spec{diss.spec.id};

            np = diss.iqc_rob.np;
            nq = diss.iqc_rob.nq;

            M_quad_rob = quad_objective_decomp(diss.iqc_rob.M, 1:np, np + (1:nq));
            [M_quad_spec, objective] = diss.spec.supply_quad(vars_spec);

            quad = obj.merge_quad(M_quad_rob, M_quad_spec);

            
            nzp = length(diss.spec.izp);
            nz = ssize(sys_cl.D, 1) - nzp;
            nwp = length(diss.spec.iwp);
            nw = ssize(sys_cl.D, 2) - nwp;

            
            Ez = eye(nz, nz+nzp);

            %constraint for the running cost
            sys_cl_run = Ez * sys_cl;
            


            %enforce L2 stability
            %the quadratic objective
            supp_b = obj.supply_block(sys_cl_run, quad);

            lam = diss.spec.weight;

            %the storage
            stor_b = obj.storage_block(sys_cl_run, quad, (1-lam) * G, G);

            %the dynamics
            [dyn_b, U_outer, V_outer] = obj.dynamics_block(sys_cl_run, quad);
            
            %wrap it all together           
            con_M = -(stor_b + supp_b + dyn_b);

            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 
            

            
            %% now do the terminal constraints
            % 
            Mterm = diss.spec.quad_terminal(vars_spec);
            % 
            % %weighting the running cost            
            quad_term = obj.merge_quad(M_quad_rob, Mterm);

            [dyn_b_term, U_outer, V_outer] = obj.dynamics_block_null(sys_cl, quad_term);

            zero_G = zeros(ssize(G, 1));
            stor_b_term = obj.storage_block(sys_cl, quad_term, (lam) * G, zero_G);



            %drop the next state evolution  
            % nt  = ssize(quad_term.U, 1);
            % nx  = ssize(G, 1);
            
            supp_b_term = obj.supply_block(sys_cl, quad_term);


            % En = eye(ssize(supp_b_term_full, 1));
            % En(:, (end-(nt+nx-1)):(end-nt)) = [];
            % 
            % 
            % supp_b_term =  En' * supp_b_term_full * En;
            % 
            % 
            % cterm_block = blkdiag(-lam* G, zeros(nw+nwp + nt));
            % 
            % n = ssize(G, 1);
            % 
            % outer_cl_right= [[ zeros(n, nw+nwp);
            %    eye(nw+nwp)], zeros(n+nw+nwp, n+nt)];
            % 
            % outer_cl_left = [quad_term.S;
            %     zeros(n, nz+nzp); %check the sign here
            %      quad_term.T];
            % 
            % 
            % % supp_b_term = obj.supply_block(sys_cl, quad_term);
            % 
            % 
            % center_cl = [sys_cl.C, sys_cl.D];
            % 
            % 
            % % center_cl = [sys_cl.A, sys_cl.B;
            % %     sys_cl.C, sys_cl.D];
            % % 
            % % outer_cl_right= [[eye(n), zeros(n, nw);
            % %     zeros(nw, n), eye(nw)], zeros(n+nw, n+nt)];
            % % 
            % 
            % 
            % term_b = outer_cl_left * center_cl * outer_cl_right; 
            % term_b_he = term_b + term_b';


            con_p2p = -(stor_b_term + supp_b_term + dyn_b_term);   
            sp2p = ssize(con_p2p, 1);
            cons = append_lmi(cons, con_p2p - obj.config.tol.M*eye(sp2p), obj.LMILAB); 

            
            %impose positivity constraint            
            cons = obj.con_terminal(G, cons, [], diss.iqc_rob);
       end        

        

        %% recovery

        function vars_new = augment_vars(obj, vars, diss, con_M)
            %AUGMENT_VARS add new variables/terms for recovery (useful for 
            %matrix elimination)                      
            vars_new = vars;            
            if obj.elimination
                vars_new.elim = con_M;            
            end

        end

        function [Ak, Bk, Ck, Dk] = recover_K_from_elim(obj, vars_rec)
            %recover the eliminated matrices in the controller   
            %    
            %Args:
            %   vars_rec: recovered variables from solver
            %
            %Returns:
            %   Ak, Bk, Ck, Dk: controller matrices
            
            if obj.elimination
                
                % https://www.sciencedirect.com/science/article/pii/0167691194000919
                %reconstruct the eliminated controller block
                M0 = vars_rec.elim.M0;
                U = vars_rec.elim.U;
                V = vars_rec.elim.V;                                

                if obj.config.syn.elimination_type == 2
                    Knull = vars_rec.elim.null;

                    M_accum = M0;

                    ns = obj.reg.ns;
                    nu = obj.sys.nu;
                    ny = obj.sys.ny;

                    nxi = size(vars_rec.diss.GY, 1);

                    K_block = zeros(ns + nxi + nu, nxi + ny);

                    [U_coord, V_coord] = obj.get_K_tri_basis(nxi);

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

                    Ak = K_block(nu + (1:nxi), 1:nxi);
                    Bk = K_block(nu + (1:nxi), (nxi+1):end);

                    Ck1 = K_block(nu + nxi + (1:ns), 1:nxi);
                    Dk1 = K_block(nu + nxi + (1:ns), (nxi+1):end);

                    Ck2 = K_block(1:nu, 1:nxi);
                    Dk2 = K_block(1:nu, (nxi+1):end);

                    Ck = [Ck1; Ck2];
                    Dk = [Dk1; Dk2];

                else
                    U = U{1};
                    V = V{2};

                    K_rem_block = basiclmi(-M0, -U, V, 'Xmin');
                
                    if obj.config.syn.elimination_type == 1
                    %get the prior entries
                        Ck2 = vars_rec.K.C;
                        Dk2 = vars_rec.K.D;
    
                        nxi = size(Ck2, 2);
    
                        %index the block
                        Ak = K_rem_block(1:nxi, 1:nxi);
                        Bk = K_rem_block(1:nxi, (nxi+1):end);
                        Ck1 = K_rem_block((nxi+1):end, 1:nxi);
                        Dk1 = K_rem_block((nxi+1):end, (nxi+1):end);
    
                        Ck = [Ck1; Ck2];
                        Dk = [Dk1; Dk2];
                    else
                        %get the prior entries
                        Bk = vars_rec.K.B;
                        Dk = vars_rec.K.D;
    
                        %index the block
                        nxi = size(Bk, 1);
                        Ak = K_rem_block(1:nxi, :);
                        Ck = K_rem_block((nxi + 1):end, :);
                    end
                end
            else
                [Ak, Bk, Ck, Dk] = recover_K_from_elim@lmi_synthesis_interface(obj, vars_rec);
            end
        end

        function el = elimination(obj)
            %ELIMINATION is the matrix elimination lemma used?
            el = obj.config.syn.elimination;   
        end

        function [sys_cl, U_cl, V_cl] = system_closed_loop(obj, P,  vars_diss, vars_reg, vars_K, rho);
            %SYSTEM_CLOSED_LOOP closed-loop matrix after nonlinear
            %transformation
            %Args:    
            %   P: IQC-filtered generalized plant 
            %   vars_diss:   variables of the problem (dissipation)
            %   vars_reg:   variables of the problem (regulator)            
            %   vars_K:   variables of the problem (controller)    
            %   rho: linear convergence rate
            %Returns:                        
            %   sys_cl:  closed-loop system dynamics
            %   U_cl:    left outer product in elimination
            %   V_cl:    right outer product in elimination
            

            if nargin < 6
                rho = 1;
            end

            rhoi = (1/rho);
            rhoiS = rhoi;
            if obj.config.gen.same_rho                
                rhoiP = 1;                
            else
                rhoiP = rhoi;
            end

            if obj.elimination
                %knock out the terms

                %get the variables
                GX = vars_diss.GX;
                GY = vars_diss.GY;

                [A, B, C, D] = ssdata(P);

                A =  rhoiP * A;
                B =  rhoiP * B;

                iu = P.index_u;
                iw = [P.index_w, P.index_wp];
               
                iy = P.index_y;
                iz = [P.index_z, P.index_zp];

                nxn = size(A, 1);                
                nu = length(iu);                
                ny = length(iy);
                nz = length(iz);
                nw = length(iw);

                
                if obj.config.syn.elimination_type == 2
                    %remove [Ak, Bk; Ck, Dk]
                    %in progress

                    ns = obj.reg.ns;
                    iu1 = iu(1:ns);         %inputs to the internal model
                    iu2 = iu((ns+1):end);   %inputs to the plant

                    nu2 = length(iu2);
                    nxi = ssize(GY, 1);

                    %closed loop without [Ak, Bk; Ck1, Dk1]

                    Acal =  [A*GY ,  A ;
                        zeros(nxi), GX*A];
                    Bcal =  [B(:, iw);
                        GX*B(:, iw)];
                    Ccal = [C(iz, :)*GY, C(iz, :) ];
                    Dcal = D(iz, iw);


                    %base outer factors
                    U_cl_base = [B(:, iu2), zeros(nxn, nxi), B(:, iu1);
                        zeros(nxi, nu2), rhoiP*eye(nxi), zeros(nxi, ns);
                        D(iz, iu2), zeros(nz, nxi), D(iz, iu1)]';
    
                    V_cl_base = [eye(nxi), zeros(nxi, nxn), zeros(nxi, nw);
                        zeros(ny, nxi), C(iy, :), D(iy, iw)];



                    %triangular decomposition
                    [U_coord, V_coord] = obj.get_K_tri_basis(nxi);

                    ntri= length(V_coord);
                    U_cl = cell(ntri+1, 1);
                    V_cl = cell(ntri+1, 1);

                    for i = 1:ntri
                        U_cl{i} = U_coord{i} * U_cl_base;
                        V_cl{i+1} = V_coord{i} * V_cl_base;
                    end



                elseif obj.config.syn.elimination_type == 1 
                    %remove [Ak, Bk; Ck1, Dk1]

                    Ck2 = vars_K.C;
                    Dk2 = vars_K.D;
                    nxi = ssize(Ck2, 2);   

                    ns = obj.reg.ns;
                    iu1 = iu(1:ns);         %inputs to the internal model
                    iu2 = iu((ns+1):end);   %inputs to the plant

                    %closed loop without [Ak, Bk; Ck1, Dk1]


                    Acal =  [A*GY + B(:, iu2)*Ck2,  A + B(:, iu2)*Dk2*C(iy, :);
                        zeros(nxi), GX*A];
                    Bcal =  [B(:, iw) + B(:, iu2)*Dk2*D(iy, iw);
                        GX*B(:, iw)];
                    Ccal = [C(iz, :)*GY+ D(iz, iu2)*Ck2, C(iz, :) + D(iz, iu2)*Dk2*C(iy, :)];
                    Dcal = D(iz, iw) + D(iz, iu2)*Dk2*D(iy, iw);


                    %outer factors
                    U_cl = [zeros(nxn, nxi), B(:, iu1);
                         rhoiP * eye(nxi), zeros(nxi, ns);
                        zeros(nz, nxi), D(iz, iu1)]';
    
                    V_cl = [eye(nxi), zeros(nxi, nxn), zeros(nxi, nw);
                        zeros(ny, nxi), C(iy, :), D(iy, iw)];


                else
                    %remove [Ak; Ck]
                    Bk = vars_K.B;                
                    Dk = vars_K.D;
                    nxi = ssize(Bk, 1);                
    
                    %closed loop without [Ak; Ck]
                    Acal = [A*GY,  A + B(:, iu)*Dk*C(iy, :);
                        zeros(nxi), GX*A + Bk*C(iy, :)];
                    Bcal = [B(:, iw) + B(:, iu)*Dk*D(iy, iw);
                        GX*B(:, iw) + Bk*D(iy, iw)];
                    Ccal = [C(iz, :)*GY, C(iz, :) + D(iz, iu)*Dk*C(iy, :)];
                    Dcal = D(iz, iw) + D(iz, iu)*Dk*D(iy, iw);                
    
    
                    
                    %outer factors
                    U_cl = [zeros(nxn, nxi), B(:, iu);
                        rhoiP * eye(nxi), zeros(nxi, nu);
                        zeros(nz, nxi), D(iz, iu)]';
    
                    V_cl = [eye(nxi), zeros(nxi, nxn+nw)];

                end

                sys_cl = sdpss(Acal, Bcal, Ccal, Dcal);
    
                
            else
                U_cl = [];
                V_cl = [];
                sys_cl = system_closed_loop@lmi_synthesis_interface(obj, P,  vars_diss, vars_reg, vars_K, rho);
            end

        end





        

    end

    
end


