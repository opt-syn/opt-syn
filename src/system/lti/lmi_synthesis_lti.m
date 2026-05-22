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

                if ~iscell(U_cl)
                    U_cl0 = U_cl;
                    V_cl0 = V_cl;
                    U_cl = {U_cl0, []};
                    V_cl = {[], V_cl0};
                end

                % if true
                    %triangular elimination

                    %from Lemma 4 of https://arxiv.org/pdf/1305.1746

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
                % else
                %     %standard elimination
                %     V_elim = V_cl * V_outer;
                %     U_elim = U_cl * U_outer;
                % 
                % 
                %     U_null = null(U_elim, 'rational');
                %     V_null = null(V_elim, 'rational');
                % 
                %     con_M_U = U_null' * con_M * U_null;
                %     con_M_V = V_null' * con_M * V_null;
                % 
                %     sMU = ssize(con_M_U,1);
                %     sMV = ssize(con_M_V,1);
                % 
                %     cons = append_lmi(cons, con_M_U - obj.config.tol.M*eye(sMU), obj.LMILAB); 
                %     cons = append_lmi(cons, con_M_V - obj.config.tol.M*eye(sMV), obj.LMILAB); 

                % end
                %store the data
                con_M_0 = con_M;

                con_M = struct;
                %the main attributes for matrix elimination
                con_M.M0 = con_M_0;
                con_M.U = U_elim;
                con_M.V = V_elim;
                con_M.null = null_accum;
                con_M.Mnull = con_M_null;
                % con_M.U = U_elim;
                % con_M.V = V_elim;

                %subsidiaries for error checking
                % con_M.U_null = U_null;
                % con_M.V_null = U_null;
                % con_M.M_U = con_M_U;
                % con_M.M_V = con_M_V;
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

        function [Ak, Bk, Ck, Dk] = recover_K_from_elim(obj, vars_rec)
            %recover the eliminated matrices in the controller            
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

                [A, B, C, D] = ssdata(P);

                
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


                    Acal = [A*GY ,  A ;
                        zeros(nxi), GX*A];
                    Bcal = [B(:, iw);
                        GX*B(:, iw)];
                    Ccal = [C(iz, :)*GY, C(iz, :) ];
                    Dcal = D(iz, iw);


                    %base outer factors
                    U_cl_base = [B(:, iu2), zeros(nxn, nxi), B(:, iu1);
                        zeros(nxi, nu2), eye(nxi), zeros(nxi, ns);
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


                    Acal = [A*GY + B(:, iu2)*Ck2,  A + B(:, iu2)*Dk2*C(iy, :);
                        zeros(nxi), GX*A];
                    Bcal = [B(:, iw) + B(:, iu2)*Dk2*D(iy, iw);
                        GX*B(:, iw)];
                    Ccal = [C(iz, :)*GY+ D(iz, iu2)*Ck2, C(iz, :) + D(iz, iu2)*Dk2*C(iy, :)];
                    Dcal = D(iz, iw) + D(iz, iu2)*Dk2*D(iy, iw);


                    %outer factors
                    U_cl = [zeros(nxn, nxi), B(:, iu1);
                        eye(nxi), zeros(nxi, ns);
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
                        eye(nxi), zeros(nxi, nu);
                        zeros(nz, nxi), D(iz, iu)]';
    
                    V_cl = [eye(nxi), zeros(nxi, nxn+nw)];

                end

                sys_cl = sdpss(Acal, Bcal, Ccal, Dcal);
    
                
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


